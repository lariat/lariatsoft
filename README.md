# LArIATSoft

This repository is meant to be set up and built on a LArIAT GPVM (lariatgpvm0X).

## Setting up the code

First log onto a LArIAT GPVM machine of your choice:

```bash
ssh -X -Y -C <USERNAME>@lariatgpvm01.fnal.gov
```

Once logged on, set up the Scientific Linux (SL7) container:

```bash
sh /exp/lariat/data/users/wforeman/SL7/start_SL7dev_jsl.sh
```

To set up a local development area for the first time:

```bash
USERDIR="/exp/lariat/app/users/<USERNAME>"
VERSION="v08_38_01"
QUALS="e17:prof"
PROJ="lariatsoft_myarea"
PROJDIR=$USERDIR/$PROJ

source /cvmfs/lariat.opensciencegrid.org/setup_lariat.sh
export MRB_PROJECT=lariatsoft

mkdir -p $PROJDIR
cd $PROJDIR
setup lariatsoft $VERSION -q $QUALS
mrb newDev
source localProducts*/setup
mrb g lariatsoft
  # To check out a specific tag or branch, uncomment these lines:
  #cd $MRB_SOURCE/lariatsoft
  #git checkout tags/$VERSION
mrbsetenv
mrb i --generator ninja
mrbslp
```

If you find that the installation stage fails, try again using only 1 CPU core:
```bash
mrb i -j 1 --generator ninja
mrbslp
```

To set up an _existing_ working area:
```bash
cd $PROJDIR
source localProducts*/setup
setup $MRB_PROJECT $MRB_PROJECT_VERSION -q $MRB_QUALS
mrbslp
```


## Setting up a specific tagged release

On the LArIAT GPVMs, can set up a specific tagged version of `lariatsoft` without having to create a development area and pull the code down from GitHub.

To see which tags are available as well as their corresponding qualifiers (ie, "e17:prof"):
```bash
ups list -aK+ lariatsoft
```

To set up a specific tag, first set up the SL7 container as described above, then do:
```bash
source /cvmfs/lariat.opensciencegrid.org/setup_lariat.sh
setup lariatsoft <tag> -q <qual>
```

Tagged release `v08_38_01a` includes all development up through February 2025.
