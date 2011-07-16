#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <getopt.h>
#include <math.h>
#include <string>
#include <map>
#include <set>
#include <iterator>
#include <iostream>

#include "datatypes_AMOS.hh"
#include "Bank_AMOS.hh"

#include "Contig_AMOS.hh"
#include "ContigEdge_AMOS.hh"
#include "Utilities_Bundler.hh"

using namespace std;
using namespace HASHMAP;
using namespace AMOS;
using namespace Bundler;

struct config {
   string bank;
   int debug;
   string overlapsName;
   string syntenyName;
};
config globals;

void printHelpText() {
   cerr << 
    "\n"
    "Filter edges indicating overlap based on an input file\n"
    "\n"
    "USAGE:\n"
    "\n"
    "FilterEdgesByOverlap -b[ank] <bank_name> -overlaps <fileName>\n"
    "The -overlaps read list of contig EID pairs that are supported by sequence overlap\n"
       << endl;
}

bool GetOptions(int argc, char ** argv) {
   int option_index = 0;      
   static struct option long_options[] = {
    {"help",               0, 0, 'h'},
    {"h",                  0, 0, 'h'},
    {"b",                  1, 0, 'b'},
    {"bank",               1, 0, 'b'},
    {"synteny",            1, 0, 's'},
    {"overlaps",           1, 0, 'o'},
    {"debug",              1, 0, 'd'},
    {0, 0, 0, 0}
  };

   globals.debug = 1;
   globals.bank = "";
   globals.overlapsName = "";
   globals.syntenyName = "";
 
   int c;
   
   while ((c = getopt_long_only(argc, argv, "", long_options, &option_index))!= -1){
      switch (c){
      case 'h':
         printHelpText();
         break;
      case 'b':
         globals.bank = string(optarg);
         break;
         break;
      case 'o':
         globals.overlapsName = string(optarg);
         break;
      case 's':
         globals.syntenyName = string(optarg);
         break;
       case 'd':
         globals.debug = atoi(optarg);
         if (globals.debug < 0) { globals.debug = 0; }
         break;
      case '?':
         return false;
      }
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

   if (globals.overlapsName == "") { // no overlaps file specified
      cerr << "A file of pair-wise overlaps must be specified " << endl;
      exit(1);
   }

   Bank_t link_bank (ContigLink_t::NCODE);
   if (! link_bank.exists(globals.bank)){
      cerr << "No edge account found in bank " << globals.bank << endl;
      exit(1);
   }
   try {
      link_bank.open(globals.bank, B_READ |B_WRITE);
   } catch (Exception_t & e) {
      cerr << "Failed to open edge account in bank " << globals.bank << ": " << endl << e << endl;
      link_bank.close();
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
   
   Bank_t contig_bank (Contig_t::NCODE);
   if (! contig_bank.exists(globals.bank)){
     cerr << "No contig account found in bank " << globals.bank << endl;
     exit(1);
   }
   try {
     contig_bank.open(globals.bank, B_READ);
   } catch (Exception_t & e) {
       cerr << "Failed to open contig account in bank " << globals.bank << ": " << endl << e << endl;
       contig_bank.close();
       exit(1);
   }

   // subset links to only those supported by true overlaps
   int counter = 0;

   ifstream overlapsFile(globals.overlapsName.c_str(), ios::in);
   if (!overlapsFile) {
      cerr << "Can't open repeats file, will not mask repeats " << globals.overlapsName << endl;
      exit(1);
   }
   string first;
   string second;
   hash_map<ID_t, set<ID_t>*, hash<ID_t>, equal_to<ID_t> > ctg2ovl;     // map from contig to edges
   while (!overlapsFile.eof()) {
     overlapsFile >> first;
     overlapsFile >> second;

     if (counter % 1000 == 0) {
        cerr << "We have read the contigs from " << counter << endl;
     }
     counter++;

     // translate to IIDs
     if (!contig_bank.existsEID(first) || !contig_bank.existsEID(second)) {
        cerr << "CONTIGS " << first << " and " << second << " dont exist" << endl;
        continue;
     }
 
     ID_t firstID = contig_bank.lookupIID(first);
     ID_t secondID = contig_bank.lookupIID(second);     
     set<ID_t>* s = ctg2ovl[firstID];
      if (s == NULL) { s = new set<ID_t>();}
      s->insert(secondID);
      ctg2ovl[firstID] = s;
      s = ctg2ovl[secondID];
      if (s == NULL) { s = new set<ID_t>();};
      s->insert(secondID);
      ctg2ovl[secondID] = s;
   }
   overlapsFile.close();
   
   counter = 0;
   hash_map<ID_t, set<ID_t>*, hash<ID_t>, equal_to<ID_t> > ctg2syn;     // map from contig to edges
   map<ID_t, bool> cte2syn;

   if (globals.syntenyName.length() > 0) {
      hash_map<ID_t, set<ID_t>*, hash<ID_t>, equal_to<ID_t> > ctg2lnk;     // map from contig to edges
      ContigEdge_t cte;
      
      for (AMOS::IDMap_t::const_iterator ci = edge_bank.getIDMap().begin(); ci; ci++) {
         edge_bank.fetch(ci->iid, cte);

         if (isBadEdge(cte)) {
            cerr << "Edge " << cte.getIID() << " ALREADY MARKED BAD, SKIPPING" << endl;
            continue;
         }
         if (cte.getIID() == 0 || cte.getContigs().first == 0 || cte.getContigs().second == 0) {
            // TODO: CTGs with links to ID 0 indicate links to/from singletons. Incorporate singletons into assembly?
           continue;
         }
         set<ID_t>* s = ctg2lnk[cte.getContigs().first];
         if (s == NULL) { s = new set<ID_t>();}
         s->insert(cte.getIID());
         ctg2lnk[cte.getContigs().first] = s;

         s = ctg2lnk[cte.getContigs().second];
         if (s == NULL) { s = new set<ID_t>();};
         s->insert(cte.getIID());
         ctg2lnk[cte.getContigs().second] = s;
      }

      ifstream syntenyFile(globals.syntenyName.c_str(), ios::in);
      if (!syntenyFile) {
         cerr << "Can't open repeats file, will not mask repeats " << globals.overlapsName << endl;
         exit(1);
      }
      Size_t size = 0;
      char adjacency;
      while (!syntenyFile.eof()) {
         syntenyFile >> first;
         syntenyFile >> second;
         syntenyFile >> size;
         syntenyFile >> adjacency;

         if (counter % 1000 == 0) {
            cerr << "We have read the contigs from " << counter << endl;
         }
         counter++;

         // translate to IIDs
         if (!contig_bank.existsEID(first) || !contig_bank.existsEID(second)) {
            cerr << "CONTIGS " << first << " and " << second << " dont exist" << endl;
            continue;
         }

         ID_t firstID = contig_bank.lookupIID(first);
         ID_t secondID = contig_bank.lookupIID(second);
         set<ID_t>* s = ctg2syn[firstID];
         if (s == NULL) { s = new set<ID_t>();}
         s->insert(secondID);
         ctg2syn[firstID] = s;
         s = ctg2syn[secondID];
         if (s == NULL) { s = new set<ID_t>();};
         s->insert(secondID);
         ctg2syn[secondID] = s;

               LinkAdjacency_t adj;
               switch(adjacency) {
                  case 'A':
                     adj = ContigEdge_t::ANTINORMAL;
                     break;
                  case 'N':
                     adj = ContigEdge_t::NORMAL;
                     break;
                  case 'I':
                     adj = ContigEdge_t::INNIE;
                     break;
                  case 'O':
                     adj = ContigEdge_t::OUTIE;
                     break;
               };

              set<ID_t> *contigLinks = ctg2lnk[firstID];
              if (contigLinks != NULL) {
               for (set<ID_t>::iterator i = contigLinks->begin(); i != contigLinks->end(); i++) {
                 edge_bank.fetch(*i, cte);
                 Size_t sizeDiff = abs(cte.getSize() - size);
                 bool storeEdge = false;

                 if (cte.getContigs().first == firstID && cte.getContigs().second == secondID && cte.getAdjacency() == adj) {
                    if ((double) sizeDiff / cte.getSD() < 5) {
                       storeEdge = true;
                    } 
                  }
                  else if (cte.getContigs().second == firstID && cte.getContigs().first == secondID) {
                     if (adj == ContigEdge_t::NORMAL) adj = ContigEdge_t::ANTINORMAL;
                     else if (adj == ContigEdge_t::ANTINORMAL) adj = ContigEdge_t::NORMAL;
                     if (cte.getAdjacency() == adj && ((double)sizeDiff / cte.getSD() < 5)) {
                        storeEdge = true;
                     }
                  } 
                  if (storeEdge == true) {
                       // add synteny link, its compatible
                       vector<ID_t> links = cte.getContigLinks();
                       // also add the synteny edge to the bank
                       ContigLink_t ctl;
                       ctl.setIID(link_bank.getIDMapSize()+1);
                       ctl.setEID(first + "_" + second);
                       ctl.setContigs(cte.getContigs());
                       ctl.setAdjacency(adj);
                       ctl.setSize(size);
                       ctl.setSD((double)size * 0.10);

                       if (link_bank.existsEID(ctl.getEID())) {
                          ID_t iid = link_bank.lookupIID(ctl.getEID());
                          link_bank.fetch(iid, ctl);
                       }
                       else {
                          link_bank.append(ctl);
                       }

                       cerr << "PROCESSING EDGE ID " << cte.getIID() << " AND ORIENT " << cte.getAdjacency() << " AND SIZE " << cte.getSize() << " KEY IS " << ctl.getEID() << " FIRST ADDING EDGE BETWEEN " << first << " AND " << second << " OF SIZE " << size << " ID " << cte.getIID() << endl;
                       links.push_back(ctl.getIID());
                       cte.setContigLinks(links);
                       edge_bank.replace(*i, cte);
                       cte2syn[cte.getIID()] = true;
                   }
                }
               }

            }
            syntenyFile.close();
      }

   ContigEdge_t cte;

   // now stream edges and mark as bad those that are not supported by overlaps
   for (AMOS::IDMap_t::const_iterator ci = edge_bank.getIDMap().begin(); ci; ci++) {
      edge_bank.fetch(ci->iid, cte);

      Contig_t ctg;
      if (cte.getSize() < 0) {
         set<ID_t>* s = ctg2ovl[cte.getContigs().first];
         string first = contig_bank.lookupEID(cte.getContigs().first);
         string second = contig_bank.lookupEID(cte.getContigs().second);

         if (s == NULL || s->find(cte.getContigs().second) == s->end()) {
            // bad edge 
            setEdgeStatus(cte, edge_bank, BAD_THRESH, false); 
cerr << "Edge " << cte.getIID() << " size: " << cte.getSize() << " FROM " << first << ", " << second << " IS BAD " << endl;
         } else {
            cerr << "Edge is supported by overlap " << cte.getIID() << " BETWEEN " << first << ", " << second << endl;
         } 
      }
   }
   flushEdgeStatus(edge_bank);

   for (AMOS::IDMap_t::const_iterator ci = edge_bank.getIDMap().begin(); ci; ci++) {
      edge_bank.fetch(ci->iid, cte);
      if (!isBadEdge(cte)) {
         string first = contig_bank.lookupEID(cte.getContigs().first);
         string second = contig_bank.lookupEID(cte.getContigs().second);
         cerr << "FINAL GOOD EDGE BETWEEN CONTIGS " << first << " " << second << " DISTANCE " << cte.getSize() << " ADJACENCY " << cte.getAdjacency();
         if (cte2syn[cte.getIID()] == true) {
            cerr << " SYNTENY SUPPORTED";
         }
         cerr << endl;
      }
   }

   for (hash_map<ID_t, set<ID_t>*, hash<ID_t>, equal_to<ID_t> >::iterator i = ctg2ovl.begin(); i != ctg2ovl.end(); i++) {
      delete (i->second);
   }

   edge_bank.close();
   link_bank.close();
   contig_bank.close();
}
