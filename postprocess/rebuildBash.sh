#!/bin/bash --login
sbatch --export=year=2015,icycle=07,nstart=1,nend=5 rebuildDomain.sh
sbatch --export=year=2015,icycle=07,nstart=6,nend=10 rebuildDomain.sh
sbatch --export=year=2015,icycle=07,nstart=11,nend=15 rebuildDomain.sh
sbatch --export=year=2015,icycle=07,nstart=16,nend=20 rebuildDomain.sh
sbatch --export=year=2015,icycle=07,nstart=21,nend=25 rebuildDomain.sh
sbatch --export=year=2015,icycle=07,nstart=26,nend=30 rebuildDomain.sh

