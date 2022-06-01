cmake_options=""

git_sources_mn=/apps/ALYA-MPIO-TOOLS/git/alya-mpio-tools.git
git_sources_github=https://github.com/dosimont/alya-mpio-tools.git
ssh_mn=mt1.bsc.es

#Reset compilation flags
CFLAGS=
CXXFLAGS=
FCFLAGS=

#Get my location
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

#Remove tmp directory

GITDIR=$DIR/git-sources

#If no sources have been found
if [ ! -d $GITDIR ]
then
  echo "No sources have been found"
  #Get the sources
  # 1. If I'm on marenostrum, try to clone from local repository:
  echo "Cloning the sources locally"
  git clone ${git_sources_mn} $GITDIR
  if [ $? -ne 0 ]
  then
    # 2. Marenostrum local repository not available, try cloning remotely:
    echo "Failed. Cloning the sources from marenostrum"
    unset SSH_AUTH_SOCK
    ssh -o PasswordAuthentication=no $USER@${ssh_mn} exit
    a=$?
    b=0
    if [ $a -eq 0 ]
    then 
      git clone $USER@${ssh_mn}:${git_sources_mn} $GITDIR
      b=$?
    fi
    if [ \( $a -ne 0 \) -o \( $b -ne 0 \) ]
    then
      # 3. Marenostrum local repository not available, try cloning from github:
      echo "Failed. Cloning the sources from github"
      git clone ${git_sources_github} $GITDIR
      if [ $? -ne 0 ]
      then
        echo "No sources have been found. Local and remote cloning unsucessful. Please download alya-mpio-tools sources manually in the following directory: $GITDIR and rerun this script"
        exit 1
      fi
    fi
  fi
fi

#Removing previous installations
rm -fr $DIR/bin $DIR/include $DIR/lib

#Going to sources
cd $GITDIR

#Updating repo
git pull

#Create and change to build directory
rm -fr build
mkdir -p build
cd build

#Configure and compile
cmake -DCMAKE_INSTALL_PREFIX=$DIR $cmake_options .. && make clean && make && make install
cd $DIR
