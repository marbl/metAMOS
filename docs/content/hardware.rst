#####################
Hardware requirements
#####################

MetAMOS was designed to work on any standard 64bit Linux
or OSX environment. To use MetAMOS for tutorial/teaching purposes, a minimum
of 16 GB RAM is recommended. To get started on real data sets a minimum of
64 GB of RAM is recommended, and up to 1 TB of RAM may be
necessary for larger datasets. In our experience, for most 50-100
million read datasets, 64-128 GB is a good place to start.

Scenario #1: running locally on large memory server 
===================================================

Suggested all-purpose build:

* 256 GB RAM  (16 x 16GB)
* 48 cores (96 HT, 4x cpu)
* 1  TB SSD temporary scratch space for running analyses
* 16 TB HDD archival space for storing analyses


Scenario #2: running on local cluster/grid
==========================================
Notes:

* Great for BLAST intense jobs
* RAM will limit the supported assemblers
* Grid job submission via SGE, others not supported
* MPI install required for Ray Meta

Scenario #3: running on cloud via Amazon EC2 High Memory 
========================================================

Recommend for best price/performance ratio -> cr1.8xlarge (Memory optimized):

* 2 X 120 GB SSD
* 32 HT x 2.8GhZ
* 244 GB RAM
* Spot instance currently at $0.361 per Hour (based on availability)
* On demand instance at $3.500 per Hour (always available)
* Reserved instance at $1.54 per Hour, $2474 upfront, approx. $2000 a month

More info on spot instances: http://aws.amazon.com/ec2/purchasing-options/spot-instances/

Or for smaller assembly jobs -> c3.8xlarge (Compute  optimized):

* 2 X 320 GB SSD
* 32 HT x 2.8GhZ
* 60 GB RAM
* On demand instance at $2.400 per Hour (spot instance same price!)

Here is a very useful cost calculator: https://www.scalyr.com/cloud/

Don't forget to account for time/cost to upload data!



