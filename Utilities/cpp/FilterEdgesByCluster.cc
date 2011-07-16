#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <getopt.h>
#include <math.h>
#include <string>
#include <map>
#include <set>
#include <queue>
#include <iterator>
#include <iostream>

#include "datatypes_AMOS.hh"
#include "Bank_AMOS.hh"
#include "BankStream_AMOS.hh"

#include "ContigIterator_AMOS.hh"

#include "Contig_AMOS.hh"
#include "ContigEdge_AMOS.hh"
#include "Scaffold_AMOS.hh"

#include "Utilities_Bundler.hh"

using namespace std;
using namespace HASHMAP;
using namespace AMOS;
using namespace Bundler;

// constants for maximum satisfied edge
static const uint32_t MAX_ITERATIONS          = 100;
static const uint32_t MAX_CLUSTERS	      = 256;
static const double_t MIN_CLUSTERED_NEIGHBORS = 0.50;
//static const double_t MIN_CLUSTERED_NEIGHBORS = 0.35;
static const double_t MIN_CLUSTER_PERCENT     = 0.75;

struct config {
   int32_t     debug;
   string      bank;
   int         maxClusterID;
   bool        removeEdges;
   HASHMAP::hash_map<ID_t, uint32_t, HASHMAP::hash<AMOS::ID_t>, HASHMAP::equal_to<AMOS::ID_t> >   clusters;  // contig clusters
};
config globals;

HASHMAP::hash_map<AMOS::ID_t, int32_t, HASHMAP::hash<AMOS::ID_t>, HASHMAP::equal_to<AMOS::ID_t> > *Bundler::cte2weight = NULL;

void printHelpText() {
   cerr << 
    "\n"
    "USAGE:\n"
    "\n"
    "FilterEdgesByCluster -b[ank] <bank_name> [-clusters fileName] <-noRemoveEdges>\n"
    "The -clusters option specifies a file containing a two-column list of contig IIDs and their cluster assignment. Cluster 0 means unassigned.\n"
    "The noRemoveEdges option will not erase edges between clusters, only output the final contig assignments.\n" 
       << endl;
}

bool GetOptions(int argc, char ** argv) {
   int option_index = 0;      
   static struct option long_options[] = {
    {"help",               0, 0, 'h'},
    {"h",                  0, 0, 'h'},
    {"b",                  1, 0, 'b'},
    {"bank",               1, 0, 'b'},
    {"noRemoveEdges",      0, 0, 'r'},
    {"clusters",           1, 0, 'c'},    
    {"debug",              1, 0, 'd'},
    {0, 0, 0, 0}
  };

   globals.debug = 1;
   globals.maxClusterID = 0;
   globals.removeEdges = true;

   int c;
   ifstream clusterFile;
   ID_t contigID;
   uint32_t clusterID;
 
   while ((c = getopt_long_only(argc, argv, "", long_options, &option_index))!= -1){
      switch (c){
      case 'h':
         printHelpText();
         break;
      case 'b':
         globals.bank = string(optarg);
         break;
      case 'c':
         // read file here
         clusterFile.open(optarg, ios::in);
         if (!clusterFile) {
            cerr << "Can't open clusters file, will not filter edges " << optarg << endl;
            break;
         }
         while (clusterFile >> contigID >> clusterID) {
           globals.clusters[contigID] = clusterID;

           if (clusterID > globals.maxClusterID) {
             globals.maxClusterID = clusterID;
           }
         }         
         clusterFile.close();
         break;
       case 'd':
         globals.debug = atoi(optarg);
         break;
      case 'r':
         globals.removeEdges = false;
         break;
      case '?':
         return false;
      }
   }
   globals.maxClusterID++;

   // print summary
   if (globals.debug >= 1) {
      cerr << "Options summary:" << endl;
      cerr << "Bank = \t" << globals.bank << endl;
   }
   
   return true;
} // GetOptions

int main(int argc, char *argv[]) {
   if (!GetOptions(argc, argv)){
      cerr << "Command line parsing failed" << endl;
      printHelpText();
      exit(1);
   }
  
   if (globals.bank == ""){ // no bank was specified
      cerr << "A bank must be specified" << endl;
      exit(1);
   }

   Bank_t contig_bank (Contig_t::NCODE);
   if (!contig_bank.exists(globals.bank)) {
      cerr << "No contig account found in bank " << globals.bank << endl;
      exit(1);
   }
   try {
      contig_bank.open(globals.bank, B_READ);
   } catch (Exception_t & e) {
      cerr << "Failed to openc ontig account in bank " << globals.bank << ": " << endl << e << endl;
      contig_bank.close();
      exit(1);
   }
 
   Bank_t edge_bank (ContigEdge_t::NCODE);
   if (! edge_bank.exists(globals.bank)){
      cerr << "No edge account found in bank " << globals.bank << endl;
      exit(1);
   }
   try {
      edge_bank.open(globals.bank, B_READ |B_WRITE);
   } catch (Exception_t & e) {
      cerr << "Failed to open edge account in bank " << globals.bank << ": " << endl << e << endl;
      edge_bank.close();
      exit(1);
   }

   ContigEdge_t cte;
   set<ID_t> toRemove;
   int maxNode = 0;

   hash_map<ID_t, set<ID_t, EdgeWeightCmp>*, hash<ID_t>, equal_to<ID_t> > ctg2lnk;     // map from contig to edges
   cte2weight = new hash_map<ID_t, int32_t, hash<ID_t>, equal_to<ID_t> >();

   hash_map<ID_t, double_t, hash<ID_t>, equal_to<ID_t> > ctg2percent;

   for (AMOS::IDMap_t::const_iterator ci = edge_bank.getIDMap().begin(); ci; ci++) {
      edge_bank.fetch(ci->iid, cte);

      if (cte.getIID() == 0 || cte.getContigs().first == 0 || cte.getContigs().second == 0) {
         // TODO: CTGs with links to ID 0 indicate links to/from singletons. Incorporate singletons into assembly?
         if (globals.debug >= 0) {
            cerr << "WARNING: link " << cte.getIID() << " (" << cte.getContigs().first << ", " << cte.getContigs().second << ") connects to a singleton, it is being ignored" << endl;
         }
         continue;
      }
      if (isBadEdge(cte)) {
         cerr << "WARNING: link " << cte.getIID() << " is marked as bad in the bank, skipping" << endl;
         continue;
      }

      (*cte2weight)[cte.getIID()] = cte.getContigLinks().size();
      set<ID_t, EdgeWeightCmp>* s = ctg2lnk[cte.getContigs().first];
      if (s == NULL) { s = new set<ID_t, EdgeWeightCmp>();}
      s->insert(cte.getIID());
      ctg2lnk[cte.getContigs().first] = s;

      s = ctg2lnk[cte.getContigs().second];
      if (s == NULL) { s = new set<ID_t, EdgeWeightCmp>();};
      s->insert(cte.getIID());
      ctg2lnk[cte.getContigs().second] = s;
cerr << "Read in edge " << cte.getContigs().first << " AND " << cte.getContigs().second << " WEIGHT " << cte.getContigLinks().size() << endl;
    }

    for (hash_map<ID_t, set<ID_t, EdgeWeightCmp>*, hash<AMOS::ID_t>, equal_to<AMOS::ID_t> >::iterator i = ctg2lnk.begin(); i != ctg2lnk.end(); i++) {
       if (i->second != NULL) {
          int total = 0;
          set<ID_t, EdgeWeightCmp>* s = i->second;
          for (set<ID_t, EdgeWeightCmp>::iterator edges = s->begin(); edges != s->end(); edges++) {
               edge_bank.fetch(*edges, cte);
               total += cte.getContigLinks().size();
          }
cerr << "NODE " << i->first << " HAS " << total << " LINKS" << endl;
       } else {
cerr << "NODE " << i->first << " HAS 0 LINKS " << endl;
       }
    }

   uint32_t iterations = 0;
   bool changed = true;
   while (changed && iterations < MAX_ITERATIONS) {
      changed = false;
      cerr << "Iteration " << iterations++ << endl;

      for (hash_map<ID_t, set<ID_t, EdgeWeightCmp>*, hash<AMOS::ID_t>, equal_to<AMOS::ID_t> >::iterator i = ctg2lnk.begin(); i != ctg2lnk.end(); i++) {
         if (i->second != NULL) {
            set<ID_t, EdgeWeightCmp>* s = i->second;
            double_t clusterPercents[MAX_CLUSTERS];
            uint32_t totalLinks = 0;
            uint32_t allLinks = 0;
            memset(clusterPercents, 0, MAX_CLUSTERS*sizeof(double_t));

            // tally up weighted counts for each cluster that our neighbors vote for
            for (set<ID_t, EdgeWeightCmp>::iterator edges = s->begin(); edges != s->end(); edges++) {
	       edge_bank.fetch(*edges, cte);
               ID_t otherID = getEdgeDestination(i->first, cte); 
               if (globals.clusters[otherID] != 0 && globals.clusters[otherID] != globals.maxClusterID) {
                  totalLinks += cte.getContigLinks().size();
                  clusterPercents[globals.clusters[otherID]] += cte.getContigLinks().size(); 
               }
               allLinks += cte.getContigLinks().size();
            }
            if (globals.debug >= 3) { cerr << "For node " << i->first << " THE CLUSTER ASSIGNMENT IS " << globals.clusters[i->first] << endl; }
            double_t max = 0;
            uint32_t maxClust = 0;
            for (uint32_t clust = 0; clust < globals.maxClusterID; clust++) {
	       clusterPercents[clust] /= totalLinks;
               if (max < clusterPercents[clust]) {
                  max = clusterPercents[clust]; 
                  maxClust = clust;
               }
            }
            if (globals.debug >= 3) { cerr << "Node " << i->first << " ASSINGED TO " << globals.clusters[i->first] << " HAS % LINKS TO CLUSTER " << ((double)totalLinks/allLinks) << " AND MAX LINKS " << max << " TO CLUSTER " << maxClust << endl; }
            ctg2percent[i->first] = (double)totalLinks/allLinks;
            if ((double)totalLinks/allLinks > MIN_CLUSTERED_NEIGHBORS) {
               if (max >= MIN_CLUSTER_PERCENT) {
                  // do nothing, it's good
               } else {
                 maxClust = globals.maxClusterID;
               }
               changed |= (maxClust != globals.clusters[i->first]);
               globals.clusters[i->first] = maxClust;
            } 
         }
      }
   }

/*
   for (hash_map<ID_t, set<ID_t, EdgeWeightCmp>*, hash<AMOS::ID_t>, equal_to<AMOS::ID_t> >::iterator i = ctg2lnk.begin(); i != ctg2lnk.end(); i++) {
      if (ctg2percent[i->first] <= MIN_CLUSTERED_NEIGHBORS) {
         globals.clusters[i->first] = globals.maxClusterID;
      }
   }
*/

   if (globals.removeEdges) {
      for (AMOS::IDMap_t::const_iterator ci = edge_bank.getIDMap().begin(); ci; ci++) {
         edge_bank.fetch(ci->iid, cte);
         if (cte.getContigs().first == 0 || cte.getContigs().second == 0) {
         } else {
	    if (globals.clusters[cte.getContigs().first] != globals.clusters[cte.getContigs().second] && globals.clusters[cte.getContigs().first] != 0 && globals.clusters[cte.getContigs().second] != 0) {
               cerr << "Removing edge between nodes " << cte.getContigs().first << " CLUSTER: " << globals.clusters[cte.getContigs().first] << " AND " << cte.getContigs().second << " CLUSTER: " << globals.clusters[cte.getContigs().second] << endl;
               toRemove.insert(cte.getIID());
            } 
            if (globals.clusters[cte.getContigs().first] == globals.maxClusterID || globals.clusters[cte.getContigs().second] == globals.maxClusterID) {
               cerr << "Removing edge between ambiguous nodes " << cte.getContigs().first << " CLUSTER: " << globals.clusters[cte.getContigs().first] << " AND " << cte.getContigs().second << " CLUSTER: " << globals.clusters[cte.getContigs().second] << endl;
               toRemove.insert(cte.getIID());
            }
         }
      }
      for (set<ID_t>::iterator i = toRemove.begin(); i != toRemove.end(); i++) {
         edge_bank.remove(*i);
      }
   }

   edge_bank.clean();
   edge_bank.close();

   // finally output updated status for all nodes
   for (AMOS::IDMap_t::const_iterator ci = contig_bank.getIDMap().begin(); ci; ci++) {
      cout << ci->iid << " " << globals.clusters[ci->iid] << endl; 
   }

   contig_bank.close();
   
   return 0;
}
