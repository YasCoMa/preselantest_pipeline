import os
import math
import time
import json

class Prediction:
    def __init__(self, path_config):
        if( os.path.isfile(path_config) ):
            with open(path_config) as g:
                self.config = json.load(g)
        else:
            print('Error: this configuration file does not exist')
    
    def _validate_pathtool(self, path_netmhcpan):
        flag=False
        if( path_netmhcpan!=None and path_netmhcpan!='' ):
            if( os.path.isdir(path_netmhcpan) ):
                self.path_tool = path_netmhcpan
                flag=True
            else:
                print('Error: this directory does not exist')
        else:
            print('Error: The path for netmhcpan')    
        
        return flag 

    def filter_proteomes(self, folder_in):
        os.system( f'rm {folder_in}*_seqs* ')
        
        for f in os.listdir(folder_in):
            name = f.split('.')[0]
            dt={}
            if(f.endswith('.faa') and f.find('_seqs.faa')==-1 ):
                g=open(folder_in+f,'r')
                for line in g:
                    l=line.replace('\n','')
                    if(l.startswith('>')):
                        k=l.replace('>','')
                        dt[k]=""
                    else:
                        dt[k]+=l
                g.close()
                 
                keys=[]
                for k in dt.keys():
                    if( len(dt[k])>=8 and len(dt[k])<=20000 ):
                        keys.append(k)
                    
                for k in keys:
                    k = k.lower().replace('/','').replace(':','').replace('+','').replace('(','').replace(')','').replace('[','').replace(']','').replace('|','_').replace(' ','_')
                    aux = list( filter( lambda x: x!='', k.split('_') ) )
                    k = '_'.join(aux)
                    with open(f'{folder_in}{name}_seqs.faa','a') as g:
                        g.write(f">{k}\n{dt[k]}\n")
                        

    def predict_epitopes(self, path_tool, folder_in, hlas=None):
        dir_out = folder_in+"results/"
        if(not os.path.isdir(dir_out)):
            os.system(f"mkdir {dir_out}")
        
        if(hlas==None):
            hlas=[ 'DRB1_0802,DRB1_0801,DRB1_0402,DRB1_0411,DRB1_0101,DRB1_0301,DRB1_0405,DRB1_0701,DRB1_0901,DRB1_1101,DRB1_1201,DRB1_1301,DRB1_1302,DRB1_1501,DRB3_0101,DRB3_0202,DRB4_0101,DRB5_0101', 'HLA-DQA10501-DQB10201,HLA-DQA10501-DQB10301,HLA-DQA10301-DQB10302,HLA-DQA10401-DQB10402,HLA-DQA10101-DQB10501,HLA-DQA10102-DQB10602', 'HLA-DPA10201-DPB10101,HLA-DPA10103-DPB10201,HLA-DPA10301-DPB10402,HLA-DPA10201-DPB10501' ]
            all_=','.join(hlas)
        else:
            all_=hlas
        
        for f in os.listdir(folder_in):
            name = f.split('.')[0]
            dt={}
            if(f.endswith('.faa') and f.find('_seqs')!=-1 ):
                name=f.split('_seqs')[0]
                print('Processing', name)
                
                path_file = f'{folder_in}{f}'
                path_out = f'{dir_out}raw_epitopes_{name}.txt'
                
                if( not os.path.isfile(path_out) ):
                    start = time.time()
                    os.system(f"{path_tool}netMHCIIpan -f {path_file} -a {all_} -BA -filter 1 -rankF 2 -rankS 1 -rankW 2  -s 1 > {path_out} ")
                    end = time.time()
                    print(name, ':', end - start)
                
    def parse_prediction_results(self, folder_in):
        dir_out = folder_in+"results/"
        if(not os.path.isdir(dir_out)):
            os.system(f"mkdir {dir_out}")
            
        gone=set()
        if( os.path.isdir( f'{dir_out}table_results.tsv') ):
            i=0
            gf=open(f'{dir_out}table_results.tsv','r')
            for line in f:
                if(i>0):
                    l=line.split('\t')
                    gone.add(l[0]+'-'+l[1])
                i+=1
            gf.close()
        else:
            gf=open(f'{dir_out}table_results.tsv','w')
            gf.write('protein\tmhc\tpeptide\tcore\tscore_ba\taffinity\trank%\n')
            gf.close()
            
        for f in os.listdir(dir_out):
            name = f.split('.')[0]
            if( f.startswith('raw_epitopes')  ):
                name=f.replace('raw_epitopes_','').split('.')[0]
                print('Processing', name)
                
                g=open(dir_out+f, 'r')
                for line in g:
                    line=line.replace('\n','')
                    if(line.find('<=WB')!=-1):
                        l = list( filter( lambda x: x!='', line.split(' ') ))
                        org=name
                        prot=l[6]
                        mhc=l[1]
                        pep=l[2]
                        core=l[4]
                        scoreba=l[10]
                        aff=l[11]
                        rank=l[12]
                        
                        if(not prot+'-'+mhc in gone):
                            with open(f'{dir_out}table_results.tsv','a') as gf:
                                gf.write(f'{prot}\t{mhc}\t{pep}\t{core}\t{scoreba}\t{aff}\t{rank}\n')
                        
                g.close()

    def parse_bepipred3_results(self, ind, inf):
        if(ind.endswith('/')):
            ind=ind[:-1]
        
        if(not os.path.isdir(f'{ind}/results')):
            os.system(f'mkdir {ind}/results')
        
        g=open(f'{ind}/results/table_results.tsv', 'w')
        g.write('protein\tpeptide\n')
        g.close()
        g=open(f'{ind}/carb_seqs.faa', 'w')
        g.close()
        
        f=open(f'{ind}/{inf}','r')
        for line in f:
            l = line.replace('\n','')
            if(l.startswith('>')):
                k = l.replace('>','').lower().replace('/','').replace(':','').replace('=','').replace('+','').replace('(','').replace(')','').replace('[','').replace(']','').replace('|','_').replace(' ','_')
                aux = list( filter( lambda x: x!='', k.split('_') ) )
                id_ = aux[0]+'_'+aux[1]
            else:
                seq=l.upper()
                with open(f'{ind}/carb_seqs.faa', 'a') as g:
                    g.write(f'>{id_}\n{seq}\n')
                    
                peps = set()
                temp=''
                for s in l:
                    if(s.isupper()):
                        temp+=s
                    else:
                        if( len(temp) >= 10 ):
                            peps.add(temp)
                        temp=''
                        
                if( len(temp) >= 10 ):
                    peps.add(temp)
                    temp=''
                    
                for p in peps:
                    with open(f'{ind}/results/table_results.tsv', 'a') as g:
                        g.write(f'{id_}\t{p}\n')
                
    def run(self, path_netmhcpan, step):
        if( path_netmhcpan!=None and self.config!=None):
            #step=self.config['step']
            for c in self.config['experiments']:
                folder_in = c['folder_in']
                
                if( c['cell_type'].lower() =='t' ):
                        if(step==0 or step==1):
                            if( self._validate_pathtool(path_netmhcpan) ):
                                hlas = c['hlas']
                                self.path_tool = path_netmhcpan
                                self.filter_proteomes(limit, folder_in)
                                self.predict_epitopes(self.path_tool, folder_in, hlas)
                                 
                        if(step==0 or step==2):
                            self.parse_prediction_results(folder_in)
                              
                if( c['cell_type'].lower() =='b' ):
                    if(step==0 or step==1):
                        self.parse_bepipred3_results(folder_in, c['bepipred3_output'] )
        else:
            print('Error: This module was not properly initialized')
                 
            
import sys

if( len(sys.argv)==3 ):
    path_config = sys.argv[1]
    step = int(sys.argv[2]) # 0 = all, 1 => filtering, 2 => prediction, 3 => parsing results
    path_netmhcpan = ''
    if( len(sys.argv)==4 ):
        path_netmhcpan = sys.argv[3] # /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/phd/drug_target_posdoc/vaxijen/new_application/netMHCIIpan-4.0/
    a=Prediction(path_config)
    a.run(path_netmhcpan, step)
else: 
    print('Error: incorrect number of passed parameters')
