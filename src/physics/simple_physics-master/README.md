# simple_physics

Simplified dry physics package, interfaced with DYNAMICO and LMDZ.

## Doxygen-generated documentation

https://ipsl.pages.in2p3.fr/projets/dynamico/simple_physics/namespaces.html

## First steps

### On Jean-Zay
Since DYNAMICO depends on XIOS, XIOS must be built. This takes more resources than
allowed on jean-zay and must be done on jean-zay-pp. However XIOS needs to be built only once.
Thus after this initial build (first step) you can work from jean-zay including to rebuild DYNAMICO.
Especially submitting jobs MUST be done from jean-zay, NOT from jean-zay-pp. 

#### GPU

* log onto jean-zay-pp to install and build XIOS and DYNAMICO 
   ~~~
    cd $WORK
    module load svn
    git clone https://gitlab.in2p3.fr/ipsl/projets/dynamico/simple_physics.git
    cd simple_physics/config/DYNAMICO
   ./DYNAMICO.sh install
   ./DYNAMICO.sh build JEANZAY_NVIDIA_ACC
  ~~~
* once XIOS and DYNAMICO are built (previous step), log onto jean-zay to submit a job
  ~~~
    cd $WORK/simple_physics/config/DYNAMICO/TEST_GPU
    sbatch -A wuu@gpu job_JEANZAY_ACC.sh
    tail -f rundir/gcm.log
    ncview rundir/physics_regular.nc
  ~~~
* rebuilding DYNAMICO after changing the code is fast and can be done on jean-zay : 
  ~~~
    cd simple_physics/config/DYNAMICO
   ./rebuild_DYNAMICO.sh
  ~~~

#### CPU

* log onto jean-zay-pp to install and build XIOS and DYNAMICO 
   ~~~
    cd $WORK
    module load svn
    git clone https://gitlab.in2p3.fr/ipsl/projets/dynamico/simple_physics.git
    cd simple_physics/config/DYNAMICO
   ./DYNAMICO.sh install
   ./DYNAMICO.sh build X64_JEANZAY
  ~~~
* once XIOS and DYNAMICO are built (previous step), log onto jean-zay to submit a job
  ~~~
    cd $WORK/simple_physics/config/DYNAMICO/TEST_trunk
    sbatch -A wuu@cpu job_JEANZAY.sh
    tail -f rundir/gcm.log
    ncview rundir/physics_regular.nc
  ~~~

## Call graph
![alt text](https://ipsl.pages.in2p3.fr/projets/dynamico/simple_physics/namespaceiniphyparam__mod_aac7c8179c3a64b82c7627f2502fb55a9_cgraph_org.svg "iniphyparam")
![alt text](https://ipsl.pages.in2p3.fr/projets/dynamico/simple_physics/namespacephyparam__mod_a44ca92c59c3a213d64b140994c4fd3d5_cgraph_org.svg "phyparam")

## Git and GitLab policies

### Branch/merge policy

Commits should not be applied directly onto the master branch. Instead :
* on your local machine, start a branch from the HEAD of `master`. Branches should be named according to :
  * feature branch : `feature/XXX`
  * bugfix branch : `hotfix/XXX`
* push your branch to the GitLab repository : `git push origin feature/XXX`
* using the GitLab interface, create a merge request. Make sure to choose the right branch onto which to merge (should be `master` by default)
* only fast-forward merges are accepted. If the master branch has changed and it is possible to rebase your branch without conflicts, the GitLab interface offers to do it.
* check that GitLab-CI pipelines pass successfully
* accept the merge request. Commits are squeezed to leave a simple history behind.  By default the source branch is deleted.

### GItLab-CI pipelines
(WIP) We currently run the following pipelines :
* All branches :
  * check that code compiles in Docker images
* `master` :
  * Generate and publish Doxygen documentation

