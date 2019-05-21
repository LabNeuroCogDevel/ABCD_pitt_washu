# ABCD
TODO: notes on project

# Notes to collaborators 
## Setup

```
mkdir ~/ABCD
git clone https://github.com/LabNeuroCogDevel/ABCD_pitt_washu/ ~/ABCD/src
```

## File hierarchy/ Data location
Proposed location to reference in scripts:

| what  | where             |
| ----- | ----------------- |
| data  | `$HOME/ABCD/data` |
| code  | `$HOME/ABCD/src`  |


Data is hosted in box, but should be downloaded locally to a path consistent across workstations. This could be symlink to system specific location, e.g. on rhea `ln -s /Volumes/Hera/Datasets/ABCD ~/ABCD/data`.



For syncing data, see [`rclone`](https://rclone.org/box/). 
  * Use like `rclone copy box:ABCD_DATA_PATH ~/ABCD/data`
  * becareful with `rclone sync`. It will delete data! 

## Using Git

see tutorials at the end of this nature [blog post](http://blogs.nature.com/naturejobs/2018/06/11/git-the-reproducibility-tool-scientists-love-to-hate/)

**RStudio ([tutorial](http://ohi-science.org/data-science-training/github.html#add-files-to-our-local-repo)) and MATLAB ([tutorial](https://osulp.github.io/git-advanced/05_integrate_matlab_and_git/index.html)) have built in git tools**

The most straight forward CLI usage looks like
```
git pull                   # get changes from everyone else -- stay up-to-date
vim 00_firstthing.bash     # edit or create a file
git add 00_firstthing.bash # tell git you want it to know about your chages (can do many times for many files)
git commit -m 'description of changes' # lock changes in
git push                   # send changes to github, where others can pull them
```

### Branches
To avoid conflicts, branch and merge: `git checkout my-cool-branch` and `git merge my-cool-branch master` 
