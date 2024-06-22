import os
import json
import pandas as pd
import numpy as np
import Levenshtein
import statistics as st

"""
inputs :
    - object: { 'saureus': { 'folder_in': 'saureus_seqs.faa', 'path_result_il': 't' } }
    - threshold rank
    - threshold number of alleles for promiscuity
    - threshold score similarity with iedb epitopes
    
scrap il inducer:
    https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php
    https://webs.iiitd.edu.in/raghava/il4pred/predict.php
    https://webs.iiitd.edu.in/raghava/il10pred/predict3.php
    http://metagenomics.iiserb.ac.in/IL17eScan/pred.php
"""

class Curation:
    def __init__(self, path_config):
        if( os.path.isfile(path_config) ):
            with open(path_config) as g:
                self.config = json.load(g)
        else:
            print('Error: this configuration file does not exist')
            
    def _filter_rank(self, df, cutoff):
        epi_out=set()
        
        for i in df.index:
            epi=df.loc[i, 'peptide']
            rank=df.loc[i, 'rank%']
            if( rank > cutoff ):
                epi_out.add(epi)
        
        return epi_out  

    def _check_promiscuity(self, folder_in, df, cutoff):
        dir_out=f"{folder_in}results/"
        epi_out=set()
        
        dt={}
        for i in df.index:
            epi=df.loc[i, 'peptide']
            mhc=df.loc[i, 'mhc']
            if( not epi in dt):
                dt[epi]=set()
            dt[epi].add(mhc)
        
        f=open( f"{dir_out}promiscuous_epitopes.tsv", "w")
        f.write('epitope\tnum_alleles\talleles\n')
        for e in dt:
            n = len(dt[e])
            if( n > cutoff ):
                epi_out.add(e)
                mhcs = ','.join( list(dt[e]) )
                f.write(f"{e}\t{len(dt[e])}\t{mhcs}\n")
        f.close() 
        
        return epi_out  
    
    def _check_overlapping_epis_violinet(self, folder_in, epitopes):
        dir_out=f"{folder_in}results/"
        epi_out=set()
        
        epis=epitopes
        
        gone=set()
        seq=''
        g=open( f'{dir_out}overlapping_violinet.tsv','w')
        g.write('protein_id\tprotein_sequence\tepitope\n')
        f=open("filters/protegen-bacterium.faa","r")
        for line in f:
            if(line!='\n'):
                l=line.replace('\n','').replace('>','')
                if(line.startswith('>')):
                    if(seq!=''):
                        for e in epis:
                            ide = e+'--'+id_
                            if(seq.find(e)!=-1 and not ide in gone ):
                                gone.add(ide)
                                epi_out.add(e)
                                g.write('%s\t%s\t%s\n' %(id_, seq, e) )
                    id_=line.split('|')[0]+'_'+line.split('|')[1]
                    seq=''
                else:
                    seq+=l
        f.close()
        g.close()
        
        if(seq!=''):
            for e in epis:
                ide = e+'--'+id_
                if(seq.find(e)!=-1 and not ide in gone ):
                    gone.add(ide)
                    epi_out.add(e)
                    g.write('%s\t%s\t%s\n' %(id_, seq, e) )
        
        return epi_out  
        
    def _filter_iedb_epitopes(self, folder_in, epitopes, cutoff):
        dir_out=f"{folder_in}results/"
        epi_out=set()
        
        epis=epitopes
        
        pubepis={}
        i=0
        f=open("filters/epitope_table_iedb.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(i>0):
                ide = l[0].split('/')[-1]
                epi = l[2]
                org = l[-4].replace('"','')
                pubepis[epi] = [ide, org]
            i+=1
        f.close()

        report={}
        details={}
        for epi in epis:
            info=[]

            similar=0
            for target in pubepis.keys():
                identity = Levenshtein.ratio(epi, target)
                if(identity >= cutoff):
                    epi_out.add(epi)
                    
                    similar+=1
                    inf = pubepis[target]
                    inf.append( str(identity*100) )
                    
                    info.append( '-'.join( inf ) )

            report[epi]=similar
            details[epi]=info

        sorted_list = sorted( report.items(), key=lambda kv: kv[1], reverse=True )  
        
        g = open( f"{dir_out}result_overlapping_iedb.tsv","w")
        g.write('epitope\tnumber_human_proteins\tepitope_iedb_info\n') 
        g.close()     
        for epi in sorted_list:
            with open( f"{dir_out}result_overlapping_iedb.tsv","a") as gf:
                gf.write("%s\t%s\t%s\n" %(epi[0], epi[1], ",".join(details[epi[0]]) ) )
                
        return epi_out
                
    def _check_human_homology(self, folder_in, epitopes):
        dir_out=f"{folder_in}results/"
        epi_out=set()
        
        humanprot=set()
        co=0
        seq=""
        mp={}
        f=open("filters/human_38_proteins.faa","r")
        for line in f:
            l=line.replace("\n","")
            if(l.find(">")!=-1):
                if(co>0):
                    humanprot.add(seq)
                    mp[seq]=id_

                seq=""
                id_=l.replace('>','')
            else:
                seq+=l
            co+=1
        f.close()
        
        if(seq!=""):
            humanprot.add(seq)
            mp[seq]=id_
            
        epis={}
        for l in epitopes:
            if(not l in epis):
                epis[l]=set()
        
        j=1
        for hs in humanprot:
            for e in epis:
                if(e in hs):
                #if(hs.find(e)!=-1):
                    epis[e].add( mp[hs] )
            j+=1
                    
        f=open( f"{dir_out}result_overlapping_human_genome.tsv","w")
        f.write('epitope\tnumber_human_proteins\tid_human_proteins\n') 
        for ep in epis.keys():
            if( len(epis[ep]) > 0):
                epi_out.add(ep)
                ids = ','.join( list(epis[ep]) )
                f.write("%s\t%i\t%s\n" %( ep, len(epis[ep]), ids ) )
        f.close()  
        
        return epi_out
        
    def generate_filtered_list(self, folder_in, threshold_sim_iedb, threshold_alleles, threshold_rank, cell):
        dir_out=f"{folder_in}results/"
        
        df=pd.read_csv( f"{dir_out}table_results.tsv", sep="\t")
        
        if( not os.path.isfile(f"{dir_out}selected_epitopes.fasta") ):
            epitopes = set(df['peptide'].unique())
            
            f=open( f"{dir_out}log_curation.tsv", "w")
            f.write( "step\tepitopes_removed\tpercentage_initial_data\n")
            
            out=set()
            if(cell.lower()=='t'):
                print('Filtering rank percentile BA')
                out = self._filter_rank( df, threshold_rank)
                print('\tRemoved ', len(out), '/', len(epitopes))
                f.write( f"rank percentile BA\t{len(out)}\t{len(out)/len(epitopes)}\n")
            
                print('Filtering promiscuity')
                aux = self._check_promiscuity( folder_in, df, threshold_alleles)
                out = out.union( aux )
                print('\tRemoved ', len(out), '/', len(epitopes))
                f.write( f"promiscuity\t{len(aux)}\t{len(aux)/len(epitopes)}\n")
            
            print('Filtering overlapping violinet')
            aux = self._check_overlapping_epis_violinet( folder_in, epitopes)
            out = out.union( aux )
            print('\tRemoved ', len(out), '/', len(epitopes))
            f.write( f"overlapping violinet protegen\t{len(aux)}\t{len(aux)/len(epitopes)}\n")
            
            if(cell.lower()=='t'):
                print('Filtering iedb high similarity epitopes')
                aux = self._filter_iedb_epitopes( folder_in, epitopes, threshold_sim_iedb)
                out = out.union( aux )
                print('\tRemoved ', len(out), '/', len(epitopes))
                f.write( f"iedb high similarity epitopes\t{len(aux)}\t{len(aux)/len(epitopes)}\n")
            
            print('Filtering epitopes in human proteins')
            aux = self._check_human_homology( folder_in, epitopes)
            out = out.union( aux )
            print('\tRemoved ', len(out), '/', len(epitopes))
            f.write( f"epitopes in human proteins\t{len(aux)}\t{len(aux)/len(epitopes)}\n")
            
            f.close()
            
            final = epitopes - out
            
            j=1
            f=open( f"{dir_out}selected_epitopes.fasta", "w")
            for e in final:
                f.write( f">peptide_{j}\n{e}\n")
                j+=1
            f.close()
            
            dt={}
            for f in os.listdir(folder_in):
                if(f.endswith('.faa') and f.find('_seqs.faa')!=-1 ):
                    g=open(folder_in+f,'r')
                    for line in g:
                        l=line.replace('\n','')
                        if(l.startswith('>')):
                            k=l.replace('>','')
                            k = k.lower().replace('/','').replace(':','').replace('+','').replace('(','').replace(')','').replace('[','').replace(']','').replace('|','_').replace(' ','_')
                            aux = list( filter( lambda x: x!='', k.split('_') ) )
                            k = aux[0]+'_'+aux[1]
                            dt[k]=""
                        else:
                            dt[k]+=l
                    g.close()
            
            f=open( f"{dir_out}selected_proteins.fasta", "w")  
            prots=set()
            for i in df.index:
                epi = df.loc[i, 'peptide']
                protein = df.loc[i, 'protein']
                
                if( (epi in final) and (not protein in prots) ):
                    prots.add(protein)
                    f.write( f">{protein}\n{dt[protein]}\n")
            f.close()   
        
        final={}
        f=open(f"{dir_out}selected_epitopes.fasta",'r')
        for line in f:
            l = line.replace('\n','')
            if( not l.startswith('>') ):
                final[l]=ide
            else:
                ide=l.replace('>','')
        f.close()
        
        dff = {}
        for i in df.index:
            peptide = df.loc[i, 'peptide']
            if(peptide in final):
                if(not peptide in dff):
                    dff[peptide] = { 'mhc': set(), 'protein': set() }
                protein = df.loc[i, 'protein']
                dff[peptide]['protein'].add(protein)
                
                if(cell.lower()=='t'):
                    allele = df.loc[i, 'mhc']
                    dff[peptide]['mhc'].add(allele)
        
        f=open( f"{dir_out}table_final_proteins_epitopes.tsv", "w")
        f.write("epitope_id\tepitope\tprotein\tmhc_alleles\n")
        for e in final:
            mhc = ','.join( list(dff[e]['mhc']) )
            protein = ','.join( list(dff[e]['protein']) )
            f.write( f"{final[e]}\t{e}\t{protein}\t{mhc}\n")
        f.close()
        
    def parse_vaxijen2_server_results(self, folder_in):
        dir_out=f"{folder_in}results/"
        
        ms=['epitope','protein']
        for m in ms:
            peps={}
            f=open( f'{dir_out}selected_{m}s.fasta','r')
            for line in f:
                l=line.replace('\n','')
                if(l!=""):
                    if(l.startswith('>')):
                        id_=l[1:]
                    else:
                        peps[id_]=[l]
            f.close()
            
            if( os.path.isfile(f'{dir_out}vaxijen_{m}_results.txt') ):
                f=open(f'{dir_out}vaxijen_{m}_results.txt','r')
                for line in f:
                    l=line.replace('\n','').split(' ')
                    if( len(l) > 4 ):
                        id_=l[0][1:]
                        
                        namecl='non-antigen'
                        cl='0'
                        if( float(l[6]) >= 0.7):
                            cl='1'
                            namecl='antigen'
                        
                        if( id_ in peps.keys() ):
                            peps[id_]+=[l[6], cl, namecl]
                f.close()
                
                f=open(f'{dir_out}vaxijen_{m}_prediction.tsv','w')
                f.write('pepid\tsequence\tvalue\tclass\tclass_name\n')
                for k in peps.keys():
                    f.write('%s\t%s\t%s\t%s\t%s\n' %(k, peps[k][0], peps[k][1], peps[k][2], peps[k][3] ) )
                f.close()  
            else:
                print(f'Error: {dir_out}vaxijen_{m}_results.txt not found')      
                
    def run(self, step):
        if( self.config!=None ):
            #step=self.config['step']
            for c in self.config['experiments']:
                folder_in = c['folder_in']
                cell = c['cell_type'].lower()
                
                threshold_sim_iedb = None
                threshold_alleles = None
                threshold_rank = None
                if(cell == 't'):
                    threshold_sim_iedb = c['threshold_sim_iedb'] 
                    threshold_alleles = c['threshold_alleles']
                    threshold_rank = c['threshold_rank']
                
                if(step==0 or step==1):
                    self.generate_filtered_list( folder_in, threshold_sim_iedb, threshold_alleles, threshold_rank, cell)
                if(step==0 or step==2):
                    self.parse_vaxijen2_server_results( folder_in)
        else:
            print('Error: This module was not properly initialized')
        
    
import sys

if( len(sys.argv)==3 ):
    path_config = sys.argv[1]
    step = int(sys.argv[2]) # 1 => filtering overlapping, 2 => parse vaxijen webpage content results
    a=Curation( path_config)
    a.run(step)
else: 
    print('Error: incorrect number of passed parameters')
    
