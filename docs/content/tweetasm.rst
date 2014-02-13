############
TweetAssembler
############

Introduction
===============

TweetAssembler is a twitter-based interface to an isolate genome assembly server powered by iMetAMOS. 

Limitations
===============

Before proceeding, its important to higlight a few important points:

- The TweetAssembler server is only able to assemble a couple of requests (at best) per day
- There exists limitations on the size of the input data. i.e. MiSeq ok, HiSeq not ok.
- TweetAssembler is nothing more than a tweet-based interface to an iMetAMOS webserver.
- Given the limited resources, job queue management is disabled. You will only be able to run a job if no other jobs are active.

Quick Start
===============

1) First, issue a request to follow @imetamos:

.. image:: f0.png

2) Next, contact the developers (treangen@gmail.com) to get your twitter account added to the `allowed accounts` list.

3) Once approved, compose a tweet to @imetamos using the following syntax:

.. code-block:: bash

    @imetamos [fastq_pair_1] [fastq_pair_2] [#ASSEMBLE] [id]

- Currently reads need to be in non-interleaved fastq format. Support exists but untested for SRA ids. 
- id simply needs to be a job-unique integer to avoid duplicate tweets in case you have to submit your job multiple times before it is ran. 

.. image:: f1.png

4) You should immediately receive a response tweet similar to:

.. image:: f2.png

5) Then simply wait for the confirmation tweet that the job was successful. 

.. image:: f3.png

6) Upon completion, you will be able to view the HTML report :

.. image:: f4.png

7) and download your assembly:

.. image:: f5.png

8) Suggestions & comments welcome!