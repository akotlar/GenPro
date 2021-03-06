use 5.10.0;
use strict;
use warnings;

package GenPro::DBManager;

our $VERSION = '0.001';

# ABSTRACT: Manages Database connection
# VERSION

#TODO: Better errors; Seem to get bad perf if copy error after each db call
#TODO: Allow missing database only if in $dbReadOnly mode
#TODO: Better singleton handling
#TODO: Completely redo API to explicitly access the envName and dbName
#separately; we now allow multiple dbs stored
#Where the default db is say the main db, at key 0

use Mouse 2;
with 'Seq::Role::Message';
use Carp qw/croak/;
# extends 'Seq::DBManager';

use Data::MessagePack;
use LMDB_File qw(:all);
use Types::Path::Tiny qw/AbsPath/;
use DDP;
use Hash::Merge::Simple qw/ merge /;
use Path::Tiny;
use Scalar::Util qw/looks_like_number/;

# We will maintain our own, internal database log for errors
use Cwd;
use Log::Fast;

# Most common error is "MDB_NOTFOUND" which isn't nec. bad.
$LMDB_File::croak_on_err = 0;

######### Public Attributes
# Flag for deleting tracks instead of inserting during patch* methods
has delete => ( is => 'rw', isa => 'Bool', default => 0, lazy => 1 );

has dryRun => ( is => 'rw', isa => 'Bool', default => 0, lazy => 1 );

has envs => (is => 'ro', isa => 'HashRef', init_arg => undef, default => sub {{}});

# has cursors => (is => 'ro', isa => 'HashRef');
#
# DBManager maintains its own, internal log, so that in a multi-user environment
# DBA can keep track of key errors
# TODO: Warning: cwd may fill up if left unchecked
my $internalLog = Log::Fast->new(
  {
    path => path( getcwd() )->child('dbManager-error.log')->stringify,
    pid  => $$,
  }
);


# instanceConfig variable holding our databases; this way can be used
# in environment where calling process never croaks
# {
#   database_dir => {
#     env => $somEnv, dbi => $someDbi
#   }
# }
####################### Static Properties ############################
my $instance;

# instanceConfig contains
# databaseDir => <Path::Tiny>
# readOnly => <Bool>
my %instanceConfig;

# Each process should open each environment only once.
# http://www.lmdb.tech/doc/starting.html
# my %envs;

# We are enforcing a singel transaction per environment for the moment, especially
# in light of the apparent LMDB_File restriction to this effect
# my %cursors;

# Can call as class method (DBManager->setDefaultDatabaseDir), or as instanceConfig method
# Prepares the class for consumption; should be run before the program can fork
# To ensure that all old data is cleared, if executing from a long-running process
sub initialize {
  my $data = @_ == 2 ? $_[1] : $_[0];

  if ($instance) {
    $internalLog->WARN("dbManager already initialized, clearing singleton");
    undef $instance;
    return;
  }

  if ( !$data->{databaseDir} ) {
    $internalLog->ERR("dbManager requires a databaseDir");
    croak;
  }

  $instanceConfig{databaseDir} = path( $data->{databaseDir} );
  if ( !$instanceConfig{databaseDir}->exists ) {
    $instanceConfig{databaseDir}->mkpath;
  }

  if ( $data->{readOnly} ) {
    $instanceConfig{readOnly} = 1;
  }

  shift;
  return __PACKAGE__->new(@_);
}

around 'new' => sub {
  my $orig = shift;
  my $self = shift;

  return $instance //= $self->$orig(@_);
};

sub BUILD {
  my $self = shift;

  # croak 'HELLO';

  # TODO: think about better way to initialize this class w.r.t databaseDir
  if ( !$instanceConfig{databaseDir} ) {
    $self->_errorWithCleanup("DBManager requires databaseDir");
  }

  if ( !$instanceConfig{databaseDir}->is_dir ) {
    $self->_errorWithCleanup('databaseDir not a directory');
  }
}

# Our packing function
#treat "1" as an integer, save more space
#treat .00012 as a single precision float, saving 4 bytes.
my $mp = Data::MessagePack->new()->prefer_integer()->prefer_float32();

################### DB Read, Write methods ############################
# Unsafe for $_[2] ; will be modified if an array is passed
# Read transactions are committed by default

# Deprecated
sub dbReadOne {
  #my ($self, $db,   $pos,  $skipCommit) = @_;
  #== $_[0], $_[1], $_[2], $_[3]  # Don't assign to avoid copy

  #It is possible not to find a database in $dbReadOnly mode (for ex: refSeq for a while didn't have chrM)
  #http://ideone.com/uzpdZ8
  #                      #$name, $dontCreate, $stringKeys, $dbNamef

  if ( !$_[1]->{db}->Alive ) {
    $_[1]->{db}->Txn = $_[1]->{env}->BeginTxn();

    # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
    $_[1]->{db}->Txn->AutoCommit(1);
  }

  $_[1]->{db}->Txn->get( $_[1]->{dbi}, $_[2], my $json );

  # Commit unless the user specifically asks not to
  #if(!$skipCommit) {
  $_[1]->{db}->Txn->commit() unless $_[3];

  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_NOTFOUND ) {
      $_[0]->_errorWithCleanup("dbRead LMDB error $LMDB_File::last_err");
      return 255;
    }

    $LMDB_File::last_err = 0;
  }

  return defined $json ? $mp->unpack($json) : undef;
}

##### Raw methods
sub dbReadOneRaw {
  #my ($self, $txn, $dbi, $pos) = @_;
  #           $_[1], $_[2], $_[3];

  $_[1]->get( $_[2], $_[3], my $json );

  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_NOTFOUND ) {
      $_[0]->_errorWithCleanup("dbReadRaw LMDB error $LMDB_File::last_err");
      return 255;
    }

    $LMDB_File::last_err = 0;
  }

  if(!defined $json) {
    return undef;
  }

  return $mp->unpack($json);
}

# Simplified dbRead: Takes a transaction from an already-opened database
# You manage your own commit
# Mutates the passed in array of positions
sub dbReadRaw {
  #my ($self, $txn, $dbi, $posAref) = @_;
  #           $_[1], $_[2], $_[3]

  if ( !ref $_[3] ) {
    goto &dbReadOneRaw;
  }

  # Modifies $posAref ($_[3]) to avoid extra allocation
  my $json;
  for my $pos ( @{ $_[3] } ) {
    #$txn       $dbi
    $_[1]->get( $_[2], $pos, $json );

    $pos = defined $json ? $mp->unpack($json) : undef;
  }

  #substantial to catch any errors
  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_NOTFOUND ) {
      #$self
      $_[0]->_errorWithCleanup("dbRead LMDB error after loop: $LMDB_File::last_err");
      return 255;
    }

    #reset the class error variable, to avoid crazy error reporting later
    $LMDB_File::last_err = 0;
  }

  #will return a single value if we were passed one value
  #return $posAref;
  return $_[3];
}

sub getTxn {
  #my ($self, $db, $flags);
  #    $_[0], $_[1], $_[2]
  # Should maybe issue warning here
  if($_[1]->{db}->Alive ) {
    $_[1]->{db}->Txn->commit();
  }

  $_[1]->{db}->Txn = $_[2] ? $_[1]->{env}->BeginTxn($_[2]) : $_[1]->{env}->BeginTxn();

  # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
  $_[1]->{db}->Txn->AutoCommit(1);

  return $_[1]->{db}->Txn;
}

# Write transactions are by default committed
sub dbPut {
  #my ( $self, $db, $pos, $data, $skipCommit ) = @_;
  #           $_[1],$_[2],$_[3], $_[4]

  if ( !$_[1]->{db}->Alive ) {
    $_[1]->{db}->Txn = $_[1]->{env}->BeginTxn();

    # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
    $_[1]->{db}->Txn->AutoCommit(1);
  }

  $_[1]->{db}->Txn->put( $_[1]->{dbi}, $_[2], $mp->pack($_[3]) );

  $_[1]->{db}->Txn->commit() unless $_[4];

  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_KEYEXIST ) {
      $_[0]->_errorWithCleanup("dbPut LMDB error: $LMDB_File::last_err");
      return 255;
    }

    #reset the class error variable, to avoid crazy error reporting later
    $LMDB_File::last_err = 0;
  }

  return 0;
}

sub dbPutRaw {
  #my ( $self, $txn,  $dbi,  $key,  $value     = @_;
      # $_[0]  $_[1], $_[2], $_[3], $_[4]
  $_[1]->put( $_[2], $_[3], $mp->pack($_[4]) );

  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_KEYEXIST && $LMDB_File::last_err != MDB_NOTFOUND ) {
      $_[0]->_errorWithCleanup("dbPut LMDB error: $LMDB_File::last_err");
      return 255;
    }

    #reset the class error variable, to avoid crazy error reporting later
    $LMDB_File::last_err = 0;
  }

  return 0;
}

#cursor version
# Read transactions are by default not committed
sub dbReadAll {
  my ( $self, $db, $skipCommit ) = @_;

  if ( !$db->{db}->Alive ) {
    $db->{db}->Txn = $db->{env}->BeginTxn();

    # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
    $db->{db}->Txn->AutoCommit(1);
  }

  # We store data in sequential, integer order
  # in all but the meta tables, which don't use this function
  # LMDB::Cursor::open($txn, $db->{dbi}, my $cursor);
  my $cursor = $db->{db}->Cursor;

  my ( $key, $value, @out );
  my $first = 1;
  while (1) {
    if ($first) {
      $cursor->_get( $key, $value, MDB_FIRST );
      $first = 0;
    } else {
      $cursor->_get( $key, $value, MDB_NEXT );
    }

    #because this error is generated right after the get
    #we want to capture it before the next iteration
    #hence this is not inside while( )
    if ( $LMDB_File::last_err == MDB_NOTFOUND ) {
      $LMDB_FILE::last_err = 0;
      last;
    }

    if ($LMDB_FILE::last_err) {
      $_[0]->_errorWithCleanup("dbReadAll LMDB error $LMDB_FILE::last_err");
      return 255;
    }

    push @out, $mp->unpack($value);
  }

  #  !$skipCommit
  $db->{db}->Txn->commit() unless $skipCommit;

  if ($LMDB_File::last_err) {
    if ( $LMDB_File::last_err != MDB_NOTFOUND ) {
      $_[0]->_errorWithCleanup("dbReadAll LMDB error at end: $LMDB_File::last_err");
      return 255;
    }

    #reset the class error variable, to avoid crazy error reporting later
    $LMDB_File::last_err = 0;
  }

  return @out ? \@out : undef;
}

################################################################################


#to store any records
#For instanceConfig, here we can store our feature name mappings, our type mappings
#whether or not a particular track has completed writing, etc
state $metaDbNamePart = '_meta';

#We allow people to update special "Meta" databases
#The difference here is that for each $databaseName, there is always
#only one meta database. Makes storing multiple meta documents in a single
#meta collection easy
#For example, users may want to store field name mappings, how many rows inserted
#whether building the database was a success, and more
sub dbReadMeta {
  my ( $self, $databaseName, $metaKey ) = @_;

  my $db = $self->_getDbi($databaseName . $metaDbNamePart, undef, {stringKeys => 1});

  # pass 1 to use string keys for meta properties
  return $self->dbReadOne($db, $metaKey);
}

#@param <String> $databaseName : whatever the user wishes to prefix the meta name with
#@param <String> $metaKey : this is our "position" in the meta database
# a.k.a the top-level key in that meta database, what type of meta data this is
#@param <HashRef|Scalar> $data : {someField => someValue} or a scalar value
sub dbPatchMeta {
  my ( $self, $databaseName, $metaKey, $data ) = @_;

  my $dbName = $databaseName . $metaDbNamePart;

  my $db = $self->_getDbi($dbName, undef, {stringKeys => 1});

  # If the user treats this metaKey as a scalar value, overwrite whatever was there
  if ( !ref $data ) {
    # undef : commit every transcation
    # 1 : use string keys
    $self->dbPut( $db, $metaKey, $data );
  } else {

    # Pass 1 to merge $data with whatever was kept at this metaKey
    # Pass 1 to use string keys for meta databases
    $self->dbPatchHash( $db, $metaKey, $data );
  }

  # Make sure that we update/sync the meta data asap, since this is critical
  # to db integrity
  $self->dbForceCommit( $db );
  return;
}

sub dbDropDatabase {
  my ( $self, $db ) = @_;

  #dbDelete returns nothing

  if ( !$db->{db}->Alive ) {
    $db->{db}->Txn = ${ $db->{env} }->BeginTxn();

    # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
    $db->{db}->Txn->AutoCommit(1);
  }

  # if only arg is not truthy, database is emptied rather than dropped
  $db->{db}->drop(1);
}

# We should really make this public and the return object from this
# to be required in any methods
# TODO: pass configuration hash, rather than N arguments
# TODO: allow opening of readOnly databases, not just environments?
# @param config: dontCreate <Bool>, readOnly <Bool>, maxDbs <Int>, mapSize <Int>, stringKeys <Bool>
sub _getDbi {

  # Exists and not defined, because in read only database we may discovercea
  # that some chromosomes don't have any data (example: hg38 refSeq chrM)
  # my ( $self, $name, $table, $onfig ) = @_;
  #       $_[0]  $_[1], $_[2],

  if (  $_[0]->envs->{ $_[1] } && $_[0]->envs->{ $_[1] }{ $_[2] || 0 } ) {
    return $_[0]->envs->{ $_[1] }{ $_[2] || 0 };
  }

  # say STDERR "Couldn't find $_[1]:$_[4]";
  my ( $self, $name, $table, $config ) = @_;

  $table //= 0;

  # MDB_NOTLS must be used, because over some # of databases,
  #
  my $flags;
  if ( $instanceConfig{readOnly} || $config->{readOnly} ) {
    #MDB_RDONLY can also be set per-transcation; it's just not mentioned
    #in the docs
    $flags = MDB_NOTLS | MDB_RDONLY | MDB_NORDAHEAD;
  } else {
    # We read synchronously during building, which is our only mixed workload
    # Without NOTLS LMDB will deadlock, at least on Mac OSX, when spawning
    # more than 501 databases each with 25 sub-databases
    # some related info here
    # https://github.com/bmatsuo/lmdb-go/issues/94
    $flags = MDB_NOTLS;
  }

  if ( !$self->envs->{$name}{env} ) {
    # Must always be at least 1;
    my $maxDbs = ( $config->{maxDbs} || 0 ) + 1;

    # Maximum db size
    # This * # of environments must be < 2^48 or maybe 2^47
    # else out of memory (64bit CPU != 64 bit addressable)
    my $mapSize = $config->{mapSize} || 16 * 1024 * 1024 * 1024;

    my $dbPath = $instanceConfig{databaseDir}->child($name);

    # Create the database, only if that is what is intended
    if ( !$dbPath->is_dir ) {
      # If dbReadOnly flag set, this database will NEVER be created during the
      # current execution cycle
      if ( $instanceConfig{readOnly} || $config->{readOnly} || $config->{dontCreate} ) {
        $self->_errorWithCleanup("Called _getDbi, but readOnly or dontCreate set and database folder doesn't exist");
      }

      $dbPath->mkpath;
    }

    $dbPath = $dbPath->stringify;

    my $env = LMDB::Env->new(
      $dbPath,
      {
        mapsize    => $mapSize,
        mode       => 0600,
        flags      => $flags,
        maxdbs     => $maxDbs,
        maxreaders => 32,
      }
    );

    if ( !$env ) {
      $self->_errorWithCleanup(
        "Failed to create environment $name:$table for $instanceConfig{databaseDir}
        beacuse of $LMDB_File::last_err"
      );
    }

    $self->envs->{$name}{env} = $env;
  }

  my $env = $self->envs->{$name}{env};

  my $txn = $env->BeginTxn();

  if ($LMDB_File::last_err) {
    $self->_errorWithCleanup(
      "Failed to associate transaction with environment $name:$table
      beacuse of $LMDB_File::last_err"
    );
  }

  my $dbFlags = 0;

  # Much faster random, somewhat faster sequential performance
  # Much smaller database size (4 byte keys, vs 6-10 byte keys)
  if ( !$config->{stringKeys} ) {
    # say STDERR "Opening $name, $namedDb as integer keys";
    $dbFlags = MDB_INTEGERKEY;
  }

  if ( !$config->{dontCreate} ) {
    $dbFlags = $dbFlags | MDB_CREATE;
  }

  my $DB = $txn->OpenDB( $table, $dbFlags );

  if ($LMDB_File::last_err) {
    $self->_errorWithCleanup(
      "Failed to open database $name:$table for
      $instanceConfig{databaseDir} beacuse of $LMDB_File::last_err"
    );
  }

  # ReadMode 1 gives memory pointer for perf reasons, not safe
  $DB->ReadMode(1);

  # Now db is open
  my $err = $txn->commit();

  # This DOES NOT fully prevent EINVAL when MDB_NOSYNC is used in conjunction with NO_TLS
  # and when we're opening very many databases
  # $env->sync(1);

  if ($err) {
    $self->_errorWithCleanup("Failed to commit open db tx because: $err");
  }

  # To prevents EINVAL when MDB_NOSYNC is used in conjunction with NO_TLS
  # our database should be opened with sync enabled.
  $env->set_flags(MDB_NOSYNC, 1);

  $table //= 0;
  $self->envs->{$name}{dbs}{$table} = {
    dbi    => $DB->dbi,
    db     => $DB,
    tflags => $flags,
    env    => $env,
    envName => $name,    # can be useful in keeping track of merged records
    tableName => $table, # can be useful in keeping track of merged records
  };

  # $envs{$name}{txn} = 0;

  return $self->envs->{$name}{dbs}{$table};
}

# TODO: test with Txn->env->sync instead of {env}->sync
sub dbForceCommit {
  my ( $self, $db ) = @_;

  if ( $db->{db}->Alive ) {
    $db->{db}->Txn->commit();
  }

  # Sync in case MDB_NOSYNC, MDB_MAPASYNC, or MDB_NOMETASYNC were enabled
  # I assume that if the user is forcing commit, they also want the state of the
  # db updated
  # sync(1) flag needed to ensure that disk buffer is flushed with MDB_NOSYNC, MAPASYNC
  $db->{env}->sync(1);

  return;
}

#  Let the environment go out of scope, so that LMDB_File can clean it up
sub closeEnv  {
  my ($self, $envName) = @_;

  delete $instance->{$envName};
}

# This can be called without instantiating Seq::DBManager, either as :: or -> class method
# @param <Seq::DBManager> $self (optional)
# @param <String> $envName (optional) : the name of a specific environment
sub cleanUp {
  my ($self, $msg) = @_;

  if ($msg) {
    $internalLog->ERR($msg);
  }

  if ( $LMDB_File::last_err
    && $LMDB_File::last_err != MDB_NOTFOUND
    && $LMDB_File::last_err != MDB_KEYEXIST )
  {
    _fatalError("dbCleanUp LMDB error: $LMDB_File::last_err");

    return 255;
  }

  undef $instance;
  return;


  # The below seems to cause issues; get 22 EINVAL
  # when opening thousands of parallel databases and calling this
  # $LMDB_File::last_err = 0;

  # if ( !%{$self->envs} && !%cursors ) {
  #   return 0;
  # }

  # if ( !%{$self->envs} && %cursors ) {
  #   _fatalError('dbManager expects no cursors if no environments opened');

  #   return 255;
  # }

  # foreach ( keys %{$self->envs} ) {
  #   $self->envs->{$_}{env}->set_flags(MDB_NOSYNC, 0);
  # }

  # # We track the unsafe stuff, just as a precaution
  # foreach ( keys %cursors ) {

  #   # Check defined because database may be empty (and will be stored as undef)
  #   if ( defined $cursors{$_} ) {

  #     for my $cursor ( values %{ $cursors{$_}{cursors} } ) {
  #       $cursor->[1]->close();
  #     }

  #     $cursors{$_}{txn}->commit();

  #     delete $cursors{$_};

  #     if ( $LMDB_File::last_err
  #       && $LMDB_File::last_err != MDB_NOTFOUND
  #       && $LMDB_File::last_err != MDB_KEYEXIST )
  #     {
  #       _fatalError("dbCleanUp LMDB error: $LMDB_File::last_err");

  #       return 255;
  #     }
  #   }
  # }

  # # May not be necessary, but easier to read, protect against references
  # # not being properly deleted
  # undef %cursors;

  # foreach ( keys %{$self->envs} ) {
  #   foreach my $db ( values %{ $self->envs->{$_}{dbs} } ) {
  #     # Check defined because database may be empty (and will be stored as undef)
  #     if ( defined $db->{db} && $db->{db}->Alive ) {
  #       $db->{db}->Txn->commit();

  #       if ( $LMDB_File::last_err
  #         && $LMDB_File::last_err != MDB_NOTFOUND
  #         && $LMDB_File::last_err != MDB_KEYEXIST )
  #       {
  #         _fatalError("dbCleanUp LMDB error during db cleanup: $LMDB_File::last_err");

  #         return 255;
  #       }
  #     }
  #   }

  #   if ( defined $self->envs->{$_}{env} ) {

  #     # Sync in case MDB_NOSYNC, MDB_MAPASYNC, or MDB_NOMETASYNC were enabled
  #     # sync(1) flag needed to ensure that disk buffer is flushed with MDB_NOSYNC, MAPASYNC
  #     # $self->envs->{$_}{env}->sync(1);
  #     $self->envs->{$_}{env}->Clean();
  #   }

  #   delete $self->envs->{$_};

  #   if ($LMDB_File::last_err) {
  #     _fatalError("dbCleanUp LMDB error during env cleanup: $LMDB_File::last_err");

  #     return 255;
  #   }
  # }

  # undef $instance;

  # return 0;
}

sub cleanAndWipeSingleton {
  my $self = shift;

  $self->cleanUp();

  undef $instance;
  undef %instanceConfig;
}

# Like DESTROY, but Moosier
sub DEMOLISH {
  my $self = shift;
  # $self->cleanUp();

  # say STDERR "DESTROY RAN";

  return;
}

# For now, we'll throw the error, until program is changed to expect error/success
# status from functions
sub _errorWithCleanup {
  my ($self, $msg) = @_;

  $self->cleanUp($msg);

  _fatalError($msg);

  return;
}

sub _fatalError {
  my $msg = @_ == 2 ? $_[1] : $_[0];

  $internalLog->ERR($msg);

  # Reset error message, not sure if this is the best way
  $LMDB_File::last_err = 0;

  __PACKAGE__->log( 'fatal', $msg );
  croak $msg;
}

### Not used

################################################################################
###### For performance reasons we may want to manage our own transactions ######
######################## WARNING: *UNSAFE* #####################################
# sub dbStartCursorTxn {
#   my ( $self, $envName, $namedDb ) = @_;

#   if ( $cursors{$envName} && $cursors{$envName}{cursors}{$namedDb} ) {
#     return $cursors{$envName}{cursors}{$namedDb};
#   }

#   #It is possible not to find a database in $dbReadOnly mode (for ex: refSeq for a while didn't have chrM)
#   #http://ideone.com/uzpdZ8

#   my $db = $self->_getDbi( $envName, undef, undef, $namedDb );

#   # TODO: Better error handling; since a cursor may be used to read or write
#   # in most cases a database not existing indicates we set readOnly or are  need to return an error if the database doesn't exist
#   if ( !$db ) {
#     $self->_errorWithCleanup(
#       "Couldn't open $envName database because it doesn't exist. readOnly is "
#       . ( $instanceConfig{readOnly} ? "set" : "not set" ) );
#     return 255;
#   }

#   # TODO: Investigate why a subtransaction isn't successfully made
#   # when using BeginTxn()
#   # If we create a txn and assign it to DB->Txn, from $db->{db}->Txn->env, before creating a txn here
#   # upon trying to use the parent transaction, we will get a crash (-30782 / BAD_TXN)
#   # no such issue arises the other way around; i.e creating this transaction, then having
#   # a normal DB->Txn created as a nested transaction
#   if ( $db->{db}->Alive ) {
#     $self->_errorWithCleanup(
#       "DB alive when calling dbStartCursorTxn,
#         LMDB_File allows only 1 txn per environment.
#         Commit DB->Txn before dbStartCursorTxn"
#     );

#     return 255;
#   }

#   $cursors{$envName} //= {
#     txn     => undef,
#     cursors => {},
#   };


#   if ( !$cursors{$envName}{txn} ) {
#     # Will throw errors saying "should be nested transaction" unlike env->BeginTxn();
#     # to protect against the above BAD_TXN issue
#     my $txn = LMDB::Txn->new( $db->{env}, $db->{tflags} );

#     $txn->AutoCommit(1);

#     $cursors{$envName}{txn} = $txn;
#   }

#   # # my $txn = $cursors{$envName}{txn};
#   # my $txn = $cursors{$envName}{txn};    #->SubTxn( $db->{tflags} );
#   #                                       # $txn->AutoCommit(1);

#   # my $DB = LMDB_File->new( $txn, $db->{dbi} );

#   # my $cursor = $DB->Cursor;


#   # # p $txn;
#   # # p $subTxn;
#   # # This means LMDB_File will not track our cursor, must close/delete manually
#   LMDB::Cursor::open( $cursors{$envName}{txn}, $db->{dbi}, my $cursor );

#   # # TODO: better error handling
#   if ( !$cursor ) {
#     $self->_errorWithCleanup("Couldn't open cursor for $_[1]");
#     return 255;
#   }

#   # # Unsafe, private LMDB_File method access but Cursor::open does not track cursors
#   # $LMDB::Txn::Txns{$$txn}{Cursors}{$$cursor} = 1;

#   $cursors{$envName}{cursors}{$namedDb} = [ undef, $cursor ];

#   # # We store data in sequential, integer order
#   # # in all but the meta tables, which don't use this function
#   # # LMDB::Cursor::open($txn, $db->{dbi}, my $cursor);
#   return $cursors{$envName}{cursors}{$namedDb};
# }

# Assumes user manages their own transactions
# Don't copy variables on the stack, since this may be called billions of times
# sub dbReadOneCursorUnsafe {

#   #my ($self, $cursor, $pos) = @_;
#   #$_[0]. $_[1].   $_[2]

#   #$cursor->[1]->_get($pos)
#   $_[1]->[1]->_get( $_[2], my $json, MDB_SET );

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND ) {

#       #$self->_errorWithCleanup
#       $_[0]->_errorWithCleanup("dbEndCursorTxn LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return defined $json ? $mp->unpack($json) : undef;
# }

# Don't copy variables on the stack, since this may be called billions of times
# Instead, modify the passed $posAref (arg 3)
# sub dbReadCursorUnsafe {

#   #my ($self, $cursor, $posAref) = @_;
#   #$_[0]. $_[1].   $_[2]

#   #foreach(@{$posAref})
#   foreach ( @{ $_[2] } ) {

#     #$cursor->[1]->_get($_, my $json, MDB_SET);
#     $_[1]->[1]->_get( $_, my $json, MDB_SET );

#     $_ = defined $json ? $mp->unpack($json) : undef;
#   }

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND ) {

#       #$self
#       $_[0]->_errorWithCleanup("dbEndCursorTxn LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   #return $posAref;
#   return $_[2];
# }

# # When you need performance, especially for genome-wide insertions
# # Be an adult, manage your own cursor
# # LMDB tells you : if you commit the cursor is closed, needs to be renewed
# # Don't copy variables from the stack, since this may be called billions of times
# sub dbPatchCursorUnsafe {

#   #my ( $self, $cursor, $chr, $trackKey, $pos, $newValue, $mergeFunc) = @_;
#   #    $_[0]. $_[1].   $_[2]. $_[3].  $_[4] $_[5].     $_[6]

#   #$cursor->[1]->_get($pos, my $json, MDB_SET);
#   $_[1]->[1]->_get( $_[4], my $json, MDB_SET );

#   my $existingValue = defined $json ? $mp->unpack($json) : [];

#   my $write = 1;
#   #                            [$trackKey]
#   if ( defined $existingValue->[ $_[3] ] ) {

#     # ($mergeFunc)
#     if ( $_[6] ) {

#       #( my $err, $existingValue->[ $trackKey ] ) =
#       #$mergeFunc->($chr, $pos, $existingValue->[$trackKey], $newValue);
#       ( my $err, $existingValue->[ $_[3] ] ) =
#       $_[6]->( $_[2], $_[4], $existingValue->[ $_[3] ], $_[5] );

#       if ($err) {
#         $_[0]->_errorWithCleanup("dbPatchCursor mergeFunc error: $err");
#         return 255;
#       }

#       # nothing to do; no value returned
#       # Leads to strange issues, Error 22 / EINVAL or environment closes
#       # if return without committing read transaction with NOTLS
#       if ( !defined $existingValue->[ $_[3] ] ) {
#         $write = 0;
#       } else {
#         $write = 1;
#       }
#     } else {

#       # No overwrite allowed by default
#       # just like dbPatch, but no overwrite option
#       # Overwrite is impossible when mergeFunc is defined
#       # TODO: remove overwrite from dbPatch
#       #return 0;
#       $write = 0;
#     }
#   } else {

#     #$existingValue->[$dbName]= $newValue;
#     $existingValue->[ $_[3] ] = $_[5];
#   }

#   if ($write) {
#     #_put as used here will not return errors if the cursor is inactive
#     # hence, "unsafe"
#     if ( defined $json ) {

#       #$cursor->[1]->_put($pos, $mp->pack($existingValue), MDB_CURRENT);
#       $_[1]->[1]->_put( $_[4], $mp->pack($existingValue), MDB_CURRENT );
#     } else {

#       #$cursor->[1]->_put($pos, $mp->pack($existingValue));
#       $_[1]->[1]->_put( $_[4], $mp->pack($existingValue) );
#     }
#   }

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND
#       && $LMDB_File::last_err != MDB_KEYEXIST )
#     {
#       #$self->_errorWithCleanup...
#       $_[0]->_errorWithCleanup("dbPatchCursorUnsafe LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

# sub dbPutCursorUnsafe {

#   #my ( $self, $cursor, $pos, $newValue) = @_;
#   #    $_[0],  $_[1]   ,$_[2],$_[3]

#   #_put as used here will not return errors if the cursor is inactive
#   # hence, "unsafe"
#   $_[1]->[1]->_put( $_[2], $mp->pack( $_[3] ) );

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND && $LMDB_File::last_err != MDB_KEYEXIST ) {
#       #$self->_errorWithCleanup...
#       $_[0]->_errorWithCleanup("dbPutCursorUnsafe LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

# # commit and close a self-managed cursor object
# # TODO: Don't close cursor if not needed
# sub dbEndCursorTxn {
#   my ( $self, $envName ) = @_;

#   if ( !defined $cursors{$envName} ) {
#     return 0;
#   }

#   for my $cursor ( values %{ $cursors{$envName}{cursors} } ) {
#     $cursor->[1]->close();
#     # Right now we don't have subtxn's, I don't really know how to fix it
#     # Without using the LMDB_File->new() way
#     # Not sure also what the benefits of subtxn's are...is it benefit of locks?
#     if(defined $cursor->[0]) {
#       $cursor->[0]->commit();
#     }
#   }

#   # The parent txn
#   $cursors{$envName}{txn}->commit();

#   delete $cursors{$envName};

#   # Allow two relatively innocuous errors, kill for anything else
#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND && $LMDB_File::last_err != MDB_KEYEXIST ) {
#       $self->_errorWithCleanup("dbEndCursorTxn LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

# Delete all values within a database; a necessity if we want to update a single track
# TODO: this may inflate database size, because very long-lived transaction
# maybe should allow to commit
# sub dbDeleteAll {
#   my ( $self, $chr, $dbName, $stringKeys ) = @_;

#   #It is possible not to find a database in dbReadOnly mode (for ex: refSeq for a while didn't have chrM)
#   #http://ideone.com/uzpdZ8
#   my $db = $self->_getDbi( $chr, 0, $stringKeys ) or return;

#   if ( !$db->{db}->Alive ) {
#     $db->{db}->Txn = ${ $db->{env} }->BeginTxn();

#     # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
#     $db->{db}->Txn->AutoCommit(1);
#   }

#   # We store data in sequential, integer order
#   # in all but the meta tables, which don't use this function
#   # LMDB::Cursor::open($txn, $db->{dbi}, my $cursor);
#   my $cursor = $db->{db}->Cursor;

#   my ( $key, $value, @out );
#   my $first = 1;
#   while (1) {
#     if ($first) {
#       $cursor->_get( $key, $value, MDB_FIRST );
#       $first = 0;
#     } else {
#       $cursor->_get( $key, $value, MDB_NEXT );
#     }

#     #because this error is generated right after the get
#     #we want to capture it before the next iteration
#     #hence this is not inside while( )
#     if ( $LMDB_File::last_err == MDB_NOTFOUND ) {
#       $LMDB_FILE::last_err = 0;
#       last;
#     }

#     if ($LMDB_FILE::last_err) {
#       $_[0]->_errorWithCleanup("dbReadAll LMDB error $LMDB_FILE::last_err");
#       return 255;
#     }

#     my $vals = $mp->unpack($value);

#     if ( $vals->[$dbName] ) {
#       $vals->[$dbName] = undef;

#       $cursor->_put( $key, $mp->pack($vals), MDB_CURRENT );
#     }
#   }

#   $db->{db}->Txn->commit();

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND ) {
#       $_[0]->_errorWithCleanup("dbReadAll LMDB error at end: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

# sub dbDeleteMeta {
#   my ( $self, $databaseName, $metaKey ) = @_;

#   my $dbName = $databaseName . $metaDbNamePart;

#   my $db = $self->_getDbi($dbName, undef, {stringKeys => 1});

#   #dbDelete returns nothing
#   # last argument means non-integer keys
#   $self->dbDelete( $databaseName . $metaDbNamePart, $metaKey, 1 );
#   return;
# }


# sub dbGetNumberOfEntries {
#   my ( $self, $db ) = @_;

#   #get database, but don't create it if it doesn't exist
#   my $db = $self->_getDbi( $chr, 1 );

#   return $db ? $db->{db}->Txn->env->stat->{entries} : 0;
# }

#Method to write a single position into the main databse
# Write transactions are by default committed
# Removed delete, overwrite capacities
# sub dbPatch {
#   my ( $self, $db,       $trackIndex, $pos,  $trackValue, $mergeFunc, $skipCommit) = @_;
#   #.   $_[0], $_[1]    , $_[2]  ,     $_[3], $_[4],       $_[5]      , $_[6]

#   die "NOT IMPLEMENTED YET";

#   if ( !$db->{db}->Alive ) {
#     $db->{db}->Txn = $db->{env}->BeginTxn();

#     # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
#     $db->{db}->Txn->AutoCommit(1);
#   }

#   my $txn = $db->{db}->Txn;

#   #zero-copy
#   #$db->{db}->Txn->get($db->{dbi}, $pos, my $json);
#   $txn->get( $db->{dbi}, $_[3], my $json );

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_NOTFOUND ) {

#       #$self->
#       $_[0]->_errorWithCleanup("dbPatch LMDB error $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   my $aref = defined $json ? $mp->unpack($json) : [];

#   my $write = 1;
#   #Undefined track values are allowed as universal-type "missing data" signal
#   #            $aref->[$trackIndex]
#   if ( defined $aref->[ $_[3] ] ) {

#     #if($mergeFunc) {
#     if ( $_[6] ) {

#       #$aref->[$trackIndex]) = $mergeFunc  ->($envName,$pos, $aref->[$trackIndex], $trackValue);
#       ( my $err, $aref->[ $_[3] ] ) = $_[6]->( $_[1], $_[4], $aref->[ $_[3] ], $_[5] );

#       if ($err) {

#         #$self
#         $_[0]->_errorWithCleanup("mergeFunc error: $err");
#         return 255;
#       }

#       # TODO: is this right? Should we set the val to undef?
#       # Nothing to update
#       if ( !defined $aref->[ $_[3] ] ) {
#         $write = 0;
#       } else {
#         $write = 1;
#       }
#     } else {
#       # If we return here, the environemnt just dissapears, not sure why
#       # Even though we don't commit, seems it should be ok?
#       # Maybe a NO_TLS issue?
#       # else {
#       #   return 0;
#       # }
#       $write = 0;
#     }

#     # Else no override
#   } else {

#     #$aref->[$trackIndex] = $trackValue
#     $aref->[ $_[3] ] = $_[5];
#   }

#   #if($self->dryRun) {
#   if ( $_[0]->dryRun ) {

#     #$self->
#     $_[0]->log( 'info', "DBManager dry run: would have dbPatch $_[1]\:$_[4]" );
#   } elsif ($write) {
#     #$txn->put($db->{dbi}, $pos, $mp->pack($aref));
#     $txn->put( $db->{dbi}, $_[4], $mp->pack($aref) );
#     $txn->commit() unless $_[7];
#   } else {
#     # say STDERR "ABORTING";
#     $txn->abort() unless $_[7];
#   }

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_KEYEXIST ) {

#       #$self->
#       $_[0]->_errorWithCleanup("dbPatch put or commit LMDB error $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

#Assumes that the posHref is
# {
#   position => {
#     feature_name => {
#       ...everything that belongs to feature_name
#     }
#   }
# }

# Method to write one key => value pair to the database, as a hash
# $pos can be any string, identifies a key within the kv database
# dataHref should be {someTrackName => someData} that belongs at $chr:$pos
# Currently only used for region tracks (currently only the Gene Track region track)
# sub dbPatchHash {
#   my ( $self, $db, $pos, $dataHref, $mergeFunc, $skipCommit, $overwrite ) = @_;

#   if ( ref $dataHref ne 'HASH' ) {
#     $self->_errorWithCleanup("dbPatchHash requires a 1-element hash of a hash");
#     return 255;
#   }

#   my $dbi = $db->{dbi};

#   if ( !$db->{db}->Alive ) {
#     $db->{db}->Txn = ${ $db->{env} }->BeginTxn();

#     # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
#     $db->{db}->Txn->AutoCommit(1);
#   }

#   #zero-copy read, don't modify $json
#   $db->{db}->Txn->get( $dbi, $pos, my $json );

#   if ( $LMDB_File::last_err && $LMDB_File::last_err != MDB_NOTFOUND ) {
#     $self->_errorWithCleanup("dbPatchHash LMDB error during get: $LMDB_File::last_err");
#     return 255;
#   }

#   $LMDB_File::last_err = 0;

#   my $href;
#   my $skip;
#   if ($json) {
#     $href = $mp->unpack($json);

#     my ( $trackKey, $trackValue ) = %{$dataHref};

#     if ( !defined $trackKey || ref $trackKey ) {
#       $self->_errorWithCleanup("dbPatchHash requires scalar trackKey");
#       return 255;
#     }

#     # Allows undefined trackValue

#     if ( defined $href->{$trackKey} ) {

#       # Deletion and insertion are mutually exclusive
#       if ( $self->delete ) {
#         delete $href->{$trackKey};
#       } elsif ($overwrite) {

#         # Merge with righthand hash taking precedence, https://ideone.com/SBbfYV
#         # Will overwrite any keys of the same name (why it's called $overwrite)
#         $href = merge $href, $dataHref;
#       } elsif ( defined $mergeFunc ) {
#         ( my $err, $href->{$trackKey} ) =
#         &$mergeFunc( $db->{envName}, $pos, $href->{$trackKey}, $trackValue );

#         if ($err) {
#           $self->_errorWithCleanup("dbPatchHash mergeFunc error: $err");
#           return 255;
#         }
#       } else {

#         # Nothing to do, value exists, we're not deleting, overwrite, or merging
#         $skip = 1;
#       }
#     } elsif ( $self->delete ) {

#       # We want to delete a non-existant key, skip
#       $skip = 1;
#     } else {
#       $href->{$trackKey} = $trackValue;
#     }
#   } elsif ( $self->delete ) {

#     # If we want to delete, and no data, there's nothing to do, skip
#     $skip = 1;
#   }

#   #insert href if we have that (only truthy if defined), or the data provided as arg
#   if ( !$skip ) {
#     if ( $self->dryRun ) {
#       $self->log( 'info', "DBManager dry run: would have dbPatchHash $chr\:$pos" );
#     } else {
#       $db->{db}->Txn->put( $db->{dbi}, $pos, $mp->pack( $href || $dataHref ) );
#     }
#   }

#   $db->{db}->Txn->commit() unless $skipCommit;

#   if ($LMDB_File::last_err) {
#     if ( $LMDB_File::last_err != MDB_KEYEXIST ) {
#       $self->_errorWithCleanup("dbPut LMDB error: $LMDB_File::last_err");
#       return 255;
#     }

#     #reset the class error variable, to avoid crazy error reporting later
#     $LMDB_File::last_err = 0;
#   }

#   return 0;
# }

# sub dbDelete {
#   my ( $self, $db, $pos ) = @_;

#   if ( $self->dryRun ) {
#     $self->log( 'info', "DBManager dry run: Would have dbDelete $chr\:$pos" );
#     return 0;
#   }

#   if ( !$db->{db}->Alive ) {
#     $db->{db}->Txn = ${ $db->{env} }->BeginTxn();

#     # not strictly necessary, but I am concerned about hard to trace abort bugs related to scope
#     $db->{db}->Txn->AutoCommit(1);
#   }

#   # Error with LMDB_File api, means $data is required as 3rd argument,
#   # even if it is undef
#   $db->{db}->Txn->del( $db->{dbi}, $pos, undef );

#   if ( $LMDB_File::last_err && $LMDB_File::last_err != MDB_NOTFOUND ) {
#     $self->_errorWithCleanup("dbDelete LMDB error: $LMDB_File::last_err");
#     return 255;
#   }

#   $LMDB_File::last_err = 0;

#   $db->{db}->Txn->commit();

#   if ($LMDB_File::last_err) {
#     $self->_errorWithCleanup("dbDelete commit LMDB error: $LMDB_File::last_err");
#     return 255;
#   }

#   #reset the class error variable, to avoid crazy error reporting later
#   $LMDB_File::last_err = 0;
#   return 0;
# }
1;
