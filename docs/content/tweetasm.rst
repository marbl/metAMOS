############
TweetAssembler
############

Getting Started
===============

Before proceeding, its important to higlight a few important points:

- The TweetAssembler server is only able to assemble a couple of requests (at best) per day
- There exists limitations on the size of the input data. i.e. MiSeq ok, HiSeq not ok.
- TweetAssembler is nothing more than a tweet-based interface to an iMetAMOS webserver.

1) First, issue a request to follow @imetamos:

.. image:: f0.png

2) Next, contact the developers (treangen@gmail.com) to get your twitter account added to the `allowed accounts` list.

3) Once approved, compose a tweet to @imetamos using the following syntax:

.. code-block:: bash

    @imetamos [fastq_pair_1] [fastq_pair_2] [#ASSEMBLE] [id]

.. image:: f1.png

4) You should immediately receive a response tweet similar to:

.. image:: f2.png

5) Then simply wait for the confirmation tweet that the job was successful. 

.. image:: f3.png

6) Upon completion, you will be able to view the HTML report :

.. image:: f4.png

7) and download your assembly:

.. image:: f5.png

8) done!