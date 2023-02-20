# This file should be "sourced" into your environment
# to set up the spack repository

DIR=$HOME/spack

git clone https://github.com/spack/spack $DIR || true

source $DIR/share/spack/setup-env.sh

spack compiler find

sload () {
  # Installs a software package and loads it into the environment
  spack install --fail-fast $@
  spack load --first $@
}

slocation () {
  # Returns the installation directory of a given software package
  HASH=`spack find --no-groups --loaded -l $@ | head -1 | awk '{print $1}'`
  spack location -i /$HASH
}
