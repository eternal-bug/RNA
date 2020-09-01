#!/bin/sh

# ================== explain ====================================================================
# This script come from the https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh
# It usually download fail.
# ===============================================================================================

cat 1>&2 <<'EOS'
Warning: Linuxbrew has been merged into Homebrew.
Please migrate to the following command:
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

EOS

exec /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
