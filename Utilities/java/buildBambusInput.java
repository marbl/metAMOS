import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class buildBambusInput {
	private static class MappingInfo {
		public String readID = null;
		
		// position on the read
		int start = 0;
		int end = 0;
	
		// position on the contig
		int contigStart = 0;
		int contigEnd = 0;
	}
	
	private static final int SUB_LEN = 25;
	private static final int NUM_READS_PER_CTG = 200;
	private static final int NUM_CTGS = 200000; 
	
	private static void readBowtieResults(String fileName, HashMap<String, ArrayList<MappingInfo>> map) throws Exception {
		BufferedReader bf = Utils.getFile(fileName, "bout");
		if (bf != null) {
			String line = null;
			int counter = 0;
			while ((line = bf.readLine()) != null) {
				if (counter % 1000000 == 0) {
					System.err.println("Read " + counter + " mapping records from " + fileName);
				}
				String[] splitLine = line.trim().split("\\s+");
				MappingInfo m = new MappingInfo();
				String contigID = splitLine[2];
				Boolean isFwd = null;
				if (splitLine[1].contains("-")) {
					isFwd = false;
				} else if (splitLine[1].contains("+")) {
					isFwd = true;
				}
				
				m.readID = splitLine[0].replaceAll("/", "_");
				m.start = 1;
				m.end = SUB_LEN;
				m.contigStart = Integer.parseInt(splitLine[3]);
				if (isFwd) {
                                   m.contigEnd = m.contigStart + splitLine[4].length() - 1;
                                } else {
                                   m.contigEnd = m.contigStart;
				   m.contigStart = m.contigEnd + splitLine[4].length() - 1;
			        }
	
				if (map.get(contigID) == null) {
					map.put(contigID, new ArrayList<MappingInfo>(NUM_READS_PER_CTG));
				}
				map.get(contigID).add(m);
				counter++;
			}
			bf.close();
		}
	}
	
	private static void outputContigRecord(PrintStream out, String contigID, String sequence, ArrayList<MappingInfo> reads) {		
		out.println("##" + contigID + " " + (reads == null ? 0 : reads.size()) + " 0 bases 00000000 checksum.");
		out.print(sequence);
		if (reads != null) {
			for (MappingInfo m : reads) {
				out.println("#" + m.readID + "(" + Math.min(m.contigEnd, m.contigStart) + ") " + (m.contigEnd >= m.contigStart ? "[]" : "[RC]") + " " + (m.end - m.start + 1) + " bases, 00000000 checksum. {" + " " + (m.contigEnd >= m.contigStart ? m.start + " " + m.end : m.end + " " + m.start) + "} <" + (m.contigEnd >= m.contigStart ? (m.contigStart+1) + " " + (m.contigEnd+1) : (m.contigEnd+1) + " " + (m.contigStart+1)) + ">");
			}
		}
	}
	private static boolean containsPrefix(HashSet<String> prefix, String name, String postfix) {
           boolean contains = false;

           for (String s : prefix) {
              if (name.contains(s + postfix)) {
                 contains = true;
                 break;
              }
           }
           return contains;
        }
        
        public static void main(String[] args) throws Exception {
	    String resultDir = System.getProperty("user.dir") + "/";
	
	    if (args.length < 2) {
	       System.err.println("Please provide an asm and read directory");
	       System.exit(1);
	    }
	
	    String asmDir = args[0];
	    String readDir = args[1];

            String outPrefix = null;
            File dir = new File(asmDir);
            if (!dir.isDirectory()) {
               System.err.println("Error, asm directory " + asmDir + " is not a directory");
               System.exit(1);
            }

            for (File fs : dir.listFiles()) {
                if (fs.getName().contains("contigs.fa")) {
                  outPrefix = fs.getName().replaceAll("\\.contigs.fa", ""); 
                }
            }
	
	    dir = new File(readDir);
	    if (!dir.isDirectory()) {
	       System.err.println("Error, read directory " + readDir + " is not a directory");
	       System.exit(1);
	    }
	
	    File[] files = dir.listFiles();
	    HashSet<String> prefixes = new HashSet<String>();
	
	    for (File fs : files) {
	    	if (fs.getName().contains("denovo_duplicates_marked.trimmed.1") || fs.getName().contains("denovo_duplicates_marked.trimmed.2")) { 
	    			prefixes.add(fs.getName().replaceAll("\\.fastq", "").replaceAll("\\.bz2", "").replaceAll("\\.1", "").replaceAll("\\.2", ""));
	    	}
	    }

            System.err.println("Prefixes for files I will read are " + prefixes);
	
	    PrintStream out = new PrintStream(new File(resultDir + outPrefix + ".fasta"));
	    for (File fs : files) {
	    	// first trim to 25bp
	    	if (containsPrefix(prefixes, fs.getName(), ".1") || containsPrefix(prefixes, fs.getName(), ".2")) { 	
System.err.println("Processing file " + fs.getName() + " FOR FASTA OUTPUT");

				BufferedReader bf = Utils.getFile(fs.getAbsolutePath(), "fastq");
				if (bf != null) {
		    		String line = null;
			        int counter = 0;
		            
		            while ((line = bf.readLine()) != null) {
		            	StringBuffer b = new StringBuffer();
		            	if (counter % 4 == 0) {
		            		out.println(line.replaceAll("@", ">").replaceAll("/", "_"));
		            	} else if ((counter - 1) % 4 == 0) {
		            		out.println(line.substring(0, SUB_LEN));
		            	}
		            	if (counter % 1000000 == 0) {
		            		System.err.println("Processed " + counter + " reads");
		            		out.flush();
		            	}
		            	counter++;
		            }
		            bf.close();
				}
	    	}
		}
	    out.close();
 
            // now generate the mate information
            String libName = resultDir + outPrefix + ".library";
            System.err.println("Building library file " + libName + " with command : cat " + resultDir + outPrefix + ".fasta |grep \">\"|grep -v \"_2\" | awk '{print substr($1, 1, length($1)-2)\"_1\\t\"substr($1, 1, length($1)-2)\"_2\\tillumina\"}' |sed s/\\>//g |sort -T `pwd` | uniq");
            Process p = Runtime.getRuntime().exec("echo library\tillumina\t100\t500");
            BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            out = new PrintStream(new File(libName));
            String line = null;
            while ((line = in.readLine()) != null) {
               out.println(line); 
            }
            in.close();
int counter = 0;
            String[] cmd = {"/bin/sh", "-c", "cat " + resultDir + outPrefix + ".fasta |grep \">\"| awk '{print substr($1, 1, length($1)-2)\"_1\\t\"substr($1, 1, length($1)-2)\"_2\\tillumina\"}' |sed s/\\>//g |sort |uniq"};
            p = Runtime.getRuntime().exec(cmd);
            in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = in.readLine()) != null) {
               out.println(line);
System.err.println("Outputting " + counter + " lines for library");
            }
            in.close();
            out.close(); 
            System.err.println("Library file built");

	    // run the bowtie aligner

/*
            System.err.println("Launching bowtie aligner: perl /fs/wrenhomes/sergek/Utils/get_singles.pl -reads " + readDir + " -assembly " + asmDir + " --threads 2");
	    p = Runtime.getRuntime().exec("perl /fs/wrenhomes/sergek/Utils/get_singles.pl -reads " + readDir + " -assembly " + asmDir + " --threads 2");
	    p.waitFor();
	    System.err.println("Bowtie finished");
*/
            HashMap<String, ArrayList<MappingInfo>> map = new HashMap<String, ArrayList<MappingInfo>>(NUM_CTGS);
            for (String prefix : prefixes) { 
		String first = resultDir + prefix + ".1.bout";
		String second  = resultDir + prefix + ".2.bout";
		if (!new File(first).exists()) {
			first = first + ".bz2";
			second = second + ".bz2";
			if (!new File(first).exists()) {
				System.err.println("Cannot find bowtie output, expected " + resultDir + prefix + ".1.bout[.bz2]");
				System.exit(1);
			}
		}
		readBowtieResults(first, map);
		readBowtieResults(second, map);
	    }
 
	    // finally run through all the contig files and build the TIGR .contig file
	    dir = new File(asmDir);
	    if (!dir.isDirectory()) {
	      System.err.println("Error, read directory " + asmDir + " is not a directory");
	      System.exit(1);
	    }
	    
	    File contigFasta = null;
	    for (File f: dir.listFiles()) {
	    	if (f.getName().contains("contigs.fa")) {
	    		contigFasta = f;
	    		break;
	    	}
	    }
	    
	    out = new PrintStream(new File(resultDir + outPrefix + ".contig"));
	    BufferedReader bf = new BufferedReader(new FileReader(contigFasta));
	    line = null;
	    String contigID = null;
	    StringBuffer contigSequence = null;
	    counter = 0;
	    while ((line = bf.readLine()) != null) {
	    	String[] splitLine = line.trim().split("\\s+"); 
	    	if (splitLine[0].startsWith(">")) {
	    		if (contigID != null) {
	    			if (counter % 10000 == 0) {
	    				System.err.println("Processed in " + counter + " contig records");    				
	    			}
	    			counter++;
	    			 
	    			outputContigRecord(out, contigID, contigSequence.toString(), map.get(contigID));
	    		}
	    		contigID = splitLine[0].replaceAll(">", "");
	    		contigSequence = new StringBuffer();
	    	} else {
	    		contigSequence.append(line + "\n");
	    	}
	    }
	    
	    if (contigID != null) {
	    	buildBambusInput.outputContigRecord(out, contigID, contigSequence.toString(), map.get(contigID));
	    }
	    
	    bf.close();
	    out.close();
		}
}
