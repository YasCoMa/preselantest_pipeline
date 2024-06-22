# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 22:39:05 2022

@author: yasmmin
"""

import os
import re
import json
import urllib.request
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
h = .02  # step size in the mesh

from sklearn import svm
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import AdaBoostClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.naive_bayes import BernoulliNB
from sklearn import tree
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.experimental import enable_hist_gradient_boosting
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix

from sklearn.cluster import KMeans

from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier

import statistics as st
from joblib import dump, load

class Implementation_vaxijen:
    descriptors={}
    
    def __init__(self):
        self.descriptors=self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open("e_table_descriptors.tsv","r")
        for line in f:
            l=line.replace("\n","").split(" ")
            aa=l[0]
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
    
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l.replace('\t','_').replace('.','_')
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')
        
        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )
        
        g.close()
        f.close()
    
    def _calculate_features(self, seq, mode):
        aac={}
        
        n=len(seq)
        
        l=8
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for i in range(5):
                for j in range(5):
                    if(mode=='auto'):
                        cond=(i==j)
                    if(mode=='cross'):
                        cond=(i!=j)
                        
                    if(cond):
                        s=0
                        
                        k=0
                        while k<(n-l):
                            aa=seq[k]
                            flag=False
                            if(aa in self.descriptors.keys()):
                                e=self.descriptors[aa]
                                flag=True
                            
                            aal=seq[k+l]
                            flag2=False
                            if(aal in self.descriptors.keys()):
                                el=self.descriptors[aal]
                                flag2=True
                            
                            if(flag and flag2):
                                s+= (e[j]*el[j]) / (n-l)
                            
                            k+=1
                        
                        aac["l"+str(l)].append(s)
        return aac                
                        
        
    def build_dataset_matrix_variance(self, mode, ide, method):
        fasta_pos=ide+"/"+method+"/dataset_pos.fasta"
        fasta_neg=ide+"/"+method+"/dataset_neg.fasta"
        
        classes={"pos": 1, "neg": 0}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 26):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, mode)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()             
                        
        
    def build_dataset_matrix_variance_testset(self, mode, ide, method, type_seq): # type_seq is epitope or entire protein
        fasta_test=ide+"/"+method+"/dataset_"+type_seq+"_test.fasta"
        
        classes={"test": 2}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 26):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, mode)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()
        
class Implementation_vaxijenModified:
    descriptors={}
    
    def __init__(self):
        self.descriptors=self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open("descriptors_pmc5549711.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            aa=l[0]
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
    
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')
        
        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )
        
        g.close()
        f.close()
    
    def build_fastas(self):
        i=0
        g=open("alternative_dataset/dataset_pos.fasta","w")
        f=open("alternative_dataset/positive_set.tsv","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            if(i>0):
                g.write(">%s\n%s\n" %( l[0]+"_"+l[2], l[1] ) )
            i+=1
        f.close()
        g.close()
        
        i=0
        g=open("alternative_dataset/dataset_neg.fasta","w")
        f=open("alternative_dataset/negative_set.tsv","r")
        for line in f:
            l=line.replace("\n","").split('\t')
            if(i>0):
                g.write(">%s\n%s\n" %( l[0]+"_"+l[2], l[1] ) )
            i+=1
        f.close()
        g.close()
        
    def _get_optimal_lag(self, fasta_pos, fasta_neg):
        seqs=set()
        f=open(fasta_pos, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        f=open(fasta_neg, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        n=len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1    
        if(i<4):
            i=4
        return i
        
    def _get_optimal_lag_test(self, fasta_test):
        seqs=set()
        f=open(fasta_test, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        n=len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1    
        if(i<4):
            i=4
        return i
    
    def _calculate_features(self, seq, max_lag):
        aac={}
        
        n=len(seq)
        
        l=max_lag
        
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for j in range(6):
                mean_=0
                for aa in seq:
                    el=self.descriptors[aa]
                    mean_+=el[j]
                mean_=mean_/n
                
                s=0
                i=0
                while i<(n-l):
                    aa=seq[i]
                    flag=False
                    if(aa in self.descriptors.keys()):
                        e=self.descriptors[aa]
                        va=e[j] - mean_
                        flag=True
                    
                    aal=seq[i+l]
                    flag2=False
                    if(aal in self.descriptors.keys()):
                        el=self.descriptors[aal]
                        vb=el[j] - mean_
                        flag2=True
                    
                    if(flag and flag2):
                        s+= (va * vb) / (n-l)
                    
                    i+=1
                        
                aac["l"+str(l)].append(s)
        return aac                
                        
        
    def build_dataset_matrix_variance(self, mode, ide, method):
        fasta_pos=ide+"/"+method+"/dataset_pos.fasta"
        fasta_neg=ide+"/"+method+"/dataset_neg.fasta"
        
        max_lag=self._get_optimal_lag(fasta_pos, fasta_neg)
        print("max lag", max_lag)
        
        classes={"pos": 1, "neg": 0}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 7):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, max_lag)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace(' ','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace(' ','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()        
                        
        
    def build_dataset_matrix_variance_testset(self, mode, ide, method, type_seq): # type_seq is epitope or entire protein
        fasta_test=ide+"/"+method+"/dataset_"+type_seq+"_test.fasta"
        
        max_lag=self._get_optimal_lag_test(fasta_test)
        print("max lag", max_lag)
        
        classes={"test": 2}
        
        g=open(ide+"/"+method+"/"+mode+"_dataset.tsv", "w")
        feas=[]
        for i in range(1, 26):
            feas.append("feature"+str(i))
            
        g.write("protein\tclass\tlag\tfeatures\n")
        
        for cl in classes.keys():
            fasta=eval("fasta_"+cl)
            if(fasta!=""):
                class_=classes[cl]
                
                self.preprocess_remove_newLine(fasta)
                
                id_=""
                f=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1], "r")
                for line in f:
                    l=line.replace("\n","")
                    if(not l.startswith(">")):
                        if(id_!=""):
                            data = eval("self._calculate_features(l, max_lag)")
                            for k in data.keys():
                                strar=[]
                                for e in data[k]:
                                    strar.append( "{:.5f}".format(e) )
                                    
                                features=','.join(strar)
                                g.write("%s\t%i\t%s\t%s\n" %(id_, class_, k, features) )
                    else:
                        try:
                            id_=l.replace(">","").split("|")[1].replace('\t','').replace('.','').replace('-','')
                        except:
                            try:
                                id_=l.replace(" ","").split("|")[0].replace('\t','').replace('.','').replace('-','')
                            except:
                                id_=""
                f.close()
            
        g.close()
        
class EvaluationModels:
    def _get_optimal_lag(self, fasta_pos, fasta_neg):
        seqs=set()
        f=open(fasta_pos, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        f=open(fasta_neg, "r")
        for line in f:
            l=line.replace("\n","")
            if(not l.startswith(">")):
                seqs.add(l)
        f.close()
        
        n=len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1        
        return i
    
    def dataset_training_evaluation_separated(self, ide, met):
        kernel = 1.0 * RBF(1.0)
        
        clfs={'svm-rbf': "svm.SVC(kernel='rbf')", 
              'sgd': 'SGDClassifier(loss="hinge", penalty="l2", max_iter=5)',
              'nearestCentroid': 'NearestCentroid()',
              'adaboost': 'AdaBoostClassifier()',
              #'gaussianProcess': 'GaussianProcessClassifier(kernel=kernel, random_state=0, n_jobs=-1)',
              'BernoulliNaiveBayes': 'BernoulliNB()',
              'DecisionTree': 'tree.DecisionTreeClassifier()',
              'GradientBoosting': "GradientBoostingClassifier( n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)",
              'NeuralNetwork': "MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1)",
              'randomForest': 'RandomForestClassifier(max_depth=5)'
        }
        
        # Separated dataset of lags 
        fg=open(ide+"/"+met+"/"+"separated_model_results.tsv","w")
        fg.write("mode\tlag\tclassifier\tmean_accuracy\tmean_f1\tmean_precision\tmean_recall\tstdev_accuracy\tstdev_f1\tstdev_precision\tstdev_recall\n")
        modes =['auto','cross']
        maxi_lag=8
        if(met=='method2' or met=='method3'):
            modes =['auto']
            fasta_pos=ide+"/"+met+"/dataset_pos.fasta"
            fasta_neg=ide+"/"+met+"/dataset_neg.fasta"
            maxi_lag=self._get_optimal_lag(fasta_pos, fasta_neg)
            
        if(met=='method3'):
            maxi_lag=2
            
        for m in modes:
            for lag in range(1, maxi_lag+1):
                df=pd.read_csv(ide+"/"+met+"/"+""+m+"_dataset.tsv", sep="\t")
                x = df[ df['lag']=='l'+str(lag) ]
                if(len(x)>1):
                    y=x.iloc[:,1]
                    X=[]
                    for i in range(len(x)):
                        aux=[]
                        values=x.iloc[i,3].split(",")
                        for v in values:
                            aux.append(float(v))
                        X.append(aux)
                        
                    for classifier in clfs:
                        clf = eval(clfs[classifier])
                        clf.fit(X, y)
                        dump(clf, ide+"/"+met+'/models/'+m+'_l'+str(lag)+'_'+classifier+'_model_trained.joblib')
                        
                        f1=cross_val_score(clf, X, y, scoring='f1', cv=10)
                        precision=cross_val_score(clf, X, y, scoring='precision', cv=10)
                        recall=cross_val_score(clf, X, y, scoring='recall', cv=10) 
                        accuracy=cross_val_score(clf, X, y, scoring='accuracy', cv=10)
                        
                        fg.write("%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %(m, 'l'+str(lag), classifier, st.mean(f1), st.mean(precision), st.mean(recall), st.mean(accuracy), st.stdev(f1), st.stdev(precision), st.stdev(recall), st.stdev(accuracy) ) )
                        
                        #for i in range(len(f1)):
                        #    fg.write(m+";l"+str(lag)+";"+classifier+";"+str(f1[i])+";"+str(precision[i])+";"+str(recall[i])+";"+str(accuracy[i])+"\n")
                        
                        ides=m+"_l"+str(lag)+"_"+classifier
        fg.close()
        
    def prepare_condensed_dataset(self, ides, met):
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        for m in modes:
            flag_stop=False
            f=open(ides+"/"+met+"/"+m+"_mapping_name_features.tsv", "w")
            f.close()
        
            index=1
            
            dfv=pd.read_csv(ides+"/"+met+"/"+m+"_dataset.tsv", sep="\t")
            """pieces=[]
            auxdfp=df[ (df['class']==1) ]
            auxdfn=df[ (df['class']==0) ]
            ns=len(auxdfp)
            if(len(auxdfp) > len(auxdfn)):
                ns=len(auxdfn)
            chosen=random.sample( list(auxdfp.index), ns)
            pieces.append(df.iloc[chosen, :])
            
            chosen=random.sample( list(auxdfn.index), ns)
            pieces.append(df.iloc[chosen, :])
                
            dfv=pd.concat( pieces )"""
            y=[]
            auxx={}
            aux=[]
            for i in range(len(dfv)):
                ide=dfv.iloc[i, 0]+'-'+str(dfv.iloc[i, 1])
                ide=ide.replace('\t','').replace(' ','').replace(':','')
                lag=dfv.iloc[i,2]
                
                if(not ide in auxx.keys()):
                    if(len(auxx.keys())>0):
                        flag_stop=True
                        
                    auxx[ide]=[]
                    y.append(dfv.iloc[i, 1])
                
                h=1
                aux=[]
                values=dfv.iloc[i,3].split(",")
                for v in values:
                    aux.append(str(v))
                    
                    if(not flag_stop):
                        with open(ides+"/"+met+"/"+m+"_mapping_name_features.tsv", "a") as f:
                            f.write("%s\t%s\n" %("f"+str(index), lag+"_f"+str(h)) )
                    h+=1
                    index+=1
                
                auxx[ide]+=aux
                
            f=open(ides+"/"+met+"/"+m+"_dataset_consolidated.tsv", "w")
            f.write("protein\tclass\tfeatures\n")
            for k in auxx.keys():
                info=k.split("-")
                f.write( ('\t'.join(info))+"\t"+(','.join(auxx[k]))+"\n")
            f.close()
            
        self.prepare_clean_features_dataset(ides, met)
    
    def prepare_clean_features_dataset(self, ide, met):
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        
        mapp={}
        for m in modes:
            mapp[m]={}
            f=open(ide+"/"+met+"/"+m+"_mapping_name_features.tsv", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                mapp[m][l[0]]=l[1]
            f.close()
        
        for m in modes:
            df=pd.read_csv(ide+"/"+met+"/"+""+m+"_dataset_consolidated.tsv", sep="\t")
            feas=[]
            X=[]
            for i in range(len(df)):
                aux=[]
                values=df.iloc[i,2].split(",")
                j=1
                for v in values:
                    aux.append(float(v))
                    if(i==0):
                        feas.append("f"+str(j))
                    j+=1
                
                X.append(aux)
                
            print(m, len(X[0]), len(feas))
            X=pd.DataFrame(X)
            #print(X)
            X.columns=feas
            X.to_csv(ide+"/"+met+"/"+m+"_features_selection_dataset.tsv", sep='\t')
            
    def execute_feature_selection(self, ide, met):
        modes =['auto','cross']
        if(met=='method2' or met=='method3'):
            modes =['auto']
        
        mapp={}
        for m in modes:
            mapp[m]={}
            f=open(ide+"/"+met+"/"+m+"_mapping_name_features.tsv", "r")
            for line in f:
                l=line.replace("\n","").split("\t")
                mapp[m][l[0]]=l[1]
            f.close()
        
        clfs={'randomForest': 'RandomForestClassifier(max_depth=5)',
              'GradientBoosting': "GradientBoostingClassifier( n_estimators=100, learning_rate=1.0, max_depth=5, random_state=0)"
        }
        
        X=[]
        Y=[]
        fg=open(ide+"/"+met+"/consolidated_featureSelection_model_results.tsv","w")
        fg.write("selected_features\tmode\tclassifier\tmean_accuracy\tmean_f1\tmean_precision\tmean_recall\tstdev_accuracy\tstdev_f1\tstdev_precision\tstdev_recall\n")
        
        for m in modes:
            df=pd.read_csv(ide+"/"+met+"/"+""+m+"_dataset_consolidated.tsv", sep="\t")
            pieces=[]
            auxdfp=df[ (df['class']==1) ]
            auxdfn=df[ (df['class']==0) ]
            ns=len(auxdfp)
            if(len(auxdfp) > len(auxdfn)):
                ns=len(auxdfn)
            chosen=random.sample( list(auxdfp.index), ns)
            pieces.append(df.iloc[chosen, :])
            
            chosen=random.sample( list(auxdfn.index), ns)
            pieces.append(df.iloc[chosen, :])
                
            dfv=pd.concat( pieces )
            Y=dfv.iloc[:,1]
            X=[]
            feas=[]
            for i in range(len(dfv)):
                aux=[]
                values=dfv.iloc[i,2].split(",")
                j=1
                for v in values:
                    aux.append(float(v))
                    if(i==0):
                        feas.append("f"+str(j))
                    j+=1
                
                X.append(aux)
                
            print(m, len(X[0]), len(feas))
            X=pd.DataFrame(X)
            #print(X)
            X.columns=feas
            X.to_csv(ide+"/"+met+"/"+m+"_features_selection_dataset.tsv", sep='\t')
            
            for classifier in clfs.keys():
                print("-----", m, classifier)
                
                #forest = RandomForestClassifier(n_jobs=-1,  max_depth=5)
                forest=eval(clfs[classifier])
                forest.fit(X, Y)
                feat_selector = BorutaPy(forest, n_estimators='auto', random_state=1)
                feat_selector.fit(X.to_numpy(), Y)
                
                feature_ranks = list(zip(X.columns, feat_selector.ranking_, feat_selector.support_))
                sels=[]
                f=open(ide+"/"+met+"/"+classifier+"_"+m+"_features_selection_result.tsv","w")
                # iterate through and print out the results
                for feat in feature_ranks:
                    f.write('Feature: {:<25} Rank: {},  Keep: {} \n'.format(feat[0], feat[1], feat[2]))
                    if(feat[2]):
                        sels.append(mapp[m][feat[0]])
                f.close()
                    
                X_filtered = feat_selector.transform(X.to_numpy())
                if( len(X_filtered[0,:]) > 0):
                    X_test=X_filtered[:20000,:]
                    y_test=Y[:20000]
                    clf = AdaBoostClassifier()
                    clf.fit(X_filtered, Y)
                    dump(clf, ide+"/"+met+'/models/consolidated_'+m+'_'+classifier+'_model_trained.joblib')
                    
                    f1=cross_val_score(clf, X, Y, scoring='f1', cv=10)
                    precision=cross_val_score(clf, X, Y, scoring='precision', cv=10)
                    recall=cross_val_score(clf, X, Y, scoring='recall', cv=10) 
                    accuracy=cross_val_score(clf, X, Y, scoring='accuracy', cv=10)
                    
                    fg.write("%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %((','.join(sels)), m, classifier, st.mean(f1), st.mean(precision), st.mean(recall), st.mean(accuracy), st.stdev(f1), st.stdev(precision), st.stdev(recall), st.stdev(accuracy) ) )
                    
                    #predictions = clf.predict(X_test)
                    
                    #fg.write("%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\n" %((','.join(sels)), m, classifier, accuracy_score(y_test, predictions), f1_score(y_test, predictions), precision_score(y_test, predictions), recall_score(y_test, predictions) ) )
                    
        fg.close()

class AntigenicityTest:
    def __init__(self, path_config):
        if( os.path.isfile(path_config) ):
            with open(path_config) as g:
                self.config = json.load(g)
        else:
            print('Error: this configuration file does not exist')
            
    def preprocess_remove_newLine(self, fasta):
        seq=""
        id_=""
        fasta=fasta.replace("_new","")
        g=open(fasta.split(".")[0]+"_new."+fasta.split(".")[1],"w")
        f=open(fasta,"r")
        for line in f:
            l=line.replace("\n","")
            if(l.startswith(">")):
                if(seq!=""):
                    g.write("%s\n%s\n" %(id_, seq) )
                id_=l.replace('\t','_').replace('.','_').replace(" ",'')
                seq=""
            else:
                seq+=l.replace('X','').replace('Z','')

        if(seq!=""):
            g.write("%s\n%s\n" %(id_, seq) )

        g.close()
        f.close()
        
    def _init_folders_for_new_ds(self, identificador, m):
        os.system("mkdir "+identificador+"/"+m)
        os.system("mkdir "+identificador+"/"+m+"/models")
        os.system("cp "+identificador+"/dataset_* "+identificador+"/"+m+"/")
    
    def build_exps_datasets(self, folder):
        
        #folder="mtb_epitopes"    
        datasets=[ 'gram+','gram-','all_gram','bcipep','hla']
        
        methods=['method1', 'method2']
        
        ev=EvaluationModels()
        
        met1=Implementation_vaxijen()
        met2=Implementation_vaxijenModified()
        
        for ds in datasets:
            ds = folder+'/'+ds
            for m in methods:
                if(not os.path.isdir(ds+"/"+m)):
                    self._init_folders_for_new_ds(ds, m)    
                    print("----------- Preparing dataset")
                    modes =['auto','cross']
                    instance=met1
                    if(m=='method2'):
                        modes = ['auto']
                        instance=met2
                        
                    for mo in modes:
                        print("\t", ds, ' - ', m, ' - ', mo)
                        instance.build_dataset_matrix_variance(mo, ds, m)  
                    
                    print('---- evaluating separated', ds, m)
                    ev.dataset_training_evaluation_separated(ds, m)  
                    
                    print('---- preparing condensed', ds, m)
                    ev.prepare_condensed_dataset(ds, m)
                    
                    print('---- evaluating c', ds, m)
                    ev.execute_feature_selection(ds, m)        

    def select_models(self, folder):
        #folder="mtb_epitopes"    
        
        datasets=[ 'gram+','gram-','all_gram','bcipep','hla']
        #datasets=[ 'gram+','gram-','all_gram']
        methods=['method1','method2']
        aux={}
        for ds in datasets:
            if(not os.path.isdir(folder+"/"+ds+"/best_models")):
                os.system("mkdir "+folder+"/"+ds+"/best_models")
                
            for m in methods:
                df=pd.read_csv(folder+"/"+ds+"/"+m+'/separated_model_results.tsv', sep='\t')
                sdf=df.sort_values(by=['mean_f1'])
                for i in range(3):
                    mode = sdf.iloc[i,0]
                    lag = sdf.iloc[i,1]
                    clf = sdf.iloc[i,2]
                    
                    os.system('cp '+folder+"/"+ds+"/"+m+'/models/'+mode+"_"+lag+"_"+clf+"_model_trained.joblib "+folder+"/"+ds+'/best_models/'+str(i)+"-"+m+'-'+mode+"-"+lag+"_"+clf+'-separated_model.joblib')
                
                """
                df=pd.read_csv(folder+"/"+ds+"/"+m+'/consolidated_featureSelection_model_results.tsv', sep='\t')
                sdf=df.sort_values(by=['mean_f1'])
                mode = sdf.iloc[0,1]
                cls = sdf.iloc[0,2]
                os.system('cp '+folder+"/"+ds+"/"+m+'/models/consolidated_'+mode+"_"+cls+"_model_trained.joblib "+folder+"/"+ds+'/best_models/'+m+'_consolidated_model.joblib')
                """
    
    def prepare_testeset(self, folder):
        ev=EvaluationModels()
        print('eval')
        met1=Implementation_vaxijen()
        met2=Implementation_vaxijenModified()
        
        prefix=folder.split('/')[-1].split('_')[0]
        
        #folder='mtb_epitopes'
        datasets=[ 'epitope','protein']
        methods=['method1','method2']
        modes =['auto','cross']
        for ds in datasets:
            ds=prefix+'_'+ds
            for m in methods:
                if( not os.path.isdir(folder+"/"+ds+"/"+m)):
                    os.system("mkdir "+folder+"/"+ds+"/"+m)
                    os.system("cp "+folder+"/"+ds+"/*.fasta "+folder+"/"+ds+"/"+m)
                
                    modes =['auto','cross']
                    instance=met1
                    if(m=='method2'):
                        modes = ['auto']
                        instance=met2
                    
                    for mo in modes:
                        print("\t", ds, ' - ', m, ' - ', mo)
                        instance.build_dataset_matrix_variance_testset(mo, folder+'/'+ds, m, ds.replace(prefix+'_','') )
                        
                    ev.prepare_condensed_dataset(folder+'/'+ds, m)
                
    def test_models_separated(self, folder):
        prefix=folder.split('/')[-1].split('_')[0]
        
        #folder="mtb_epitopes"    
        datasets_model={ 'epitope': [ 'bcipep','hla'] , 'protein': ['gram+', 'gram-','all_gram'] }
        #datasets=[ 'gram+','gram-','all_gram']
        methods=['method1','method2']
        params={ 'method1': { 'cross': 20, 'auto': 5 }, 'method2': { 'auto': 6 }, 'method3': { 'auto': 178 } }
        passed=set()
        aux={}
        for d in datasets_model.keys():
            for ds in datasets_model[d]:
                d=prefix+'_'+d.replace(prefix+'_','')
                
                if(not os.path.isdir(folder+"/"+ds+"/best_models/results")):
                    os.system("mkdir "+folder+"/"+ds+"/best_models/results")
            
                if(not os.path.isdir(folder+"/"+ds+"/best_models/results/"+d)):
                    os.system("mkdir "+folder+"/"+ds+"/best_models/results/"+d)
                    
                for m in methods:
                    for mo in os.listdir(folder+"/"+ds+"/best_models"):
                        if(mo.find(m)!=-1):
                            name = mo.split(".")[0].replace('-separated_model','')
                            
                            index='top-'+str(int(name.split('-')[0])+1)
                            mode=name.split('-')[2]
                            lag=name.split('-')[3].split('_')[0]
                            df=pd.read_csv(folder+"/"+d+"/"+m+"/"+mode+'_dataset.tsv',sep='\t')
                            print(ds, d, m, mode, lag)
                            filt=df[ df['lag']==lag ]
                            if(len(filt)>0):
                                
                                cols=[] # making columns for the dataframe of the filtered lag dataset
                                ncols=params[m][mode]
                                for c in range(1,ncols+1):
                                    cols.append('f'+str(c))
                                    
                                feas=filt['features'] # getting float values of the features
                                X=[]
                                for f in feas:
                                    values=f.split(",")
                                    aux=[]
                                    for v in values:
                                        aux.append(float(v))
                                    X.append(aux)
                                tempdf=pd.DataFrame(X)
                                tempdf.columns=cols
                                
                                model=load(folder+"/"+ds+"/best_models/"+mo) # loading model
                                preds = model.predict(tempdf) # predicting
                                tempdf['preds_binary']=preds
                                
                                try:
                                    preds_ = model.predict_proba(tempdf) # predicting probabilities
                                    tempdf['preds_probs']=preds_
                                except:
                                    pass    
                                tempdf.to_csv(folder+"/"+ds+"/best_models/results/"+d+"/"+name+'_dataset_result.tsv',sep='\t')
                
    def test_models_feature_selection(self, folder):
        prefix=folder.split('/')[-1].split('_')[0]
        
        #folder="mtb_epitopes"    
        datasets_model={ 'epitope': [ 'bcipep','hla'] , 'protein': ['gram+', 'gram-','all_gram'] }
        #datasets=[ 'gram+','gram-','all_gram']
        methods=['method1','method2']
        params={ 'method1': { 'cross': 20, 'auto': 5 }, 'method2': { 'auto': 6 }, 'method3': { 'auto': 178 } }
        passed=set()
        aux={}
        for d in datasets_model.keys():
            for ds in datasets_model[d]:
                d=prefix+'_'+d.replace(prefix+'_','')
                if(not os.path.isdir(folder+"/"+ds+"/featSelect_results")):
                    os.system("mkdir "+folder+"/"+ds+"/featSelect_results")
            
                if(not os.path.isdir(folder+"/"+ds+"/featSelect_results/"+d)):
                    os.system("mkdir "+folder+"/"+ds+"/featSelect_results/"+d)
                    
                for m in methods:
                    df=pd.read_csv(folder+"/"+ds+"/"+m+"/consolidated_featureSelection_model_results.tsv", sep="\t") 
                    for i in range(len(df)):
                        feats=df.iloc[i,0].split(',')
                        
                        mapp={}
                        f=open(folder+"/"+ds+"/"+m+"/"+df.iloc[i,1]+"_mapping_name_features.tsv",'r')
                        for line in f:
                            l=line.replace('\n','').split('\t')
                            mapp[l[1]]=l[0]
                        f.close()
                        
                        target=pd.read_csv(folder+'/'+d+"/"+m+"/"+df.iloc[i,1]+"_features_selection_dataset.tsv", sep='\t')
                        cols=[]
                        for c in feats:
                            if(c in mapp.keys()):
                                cols.append(mapp[c])
                        
                        intersec=set(target.columns).intersection(set(cols))
                        
                        selection=list(intersec) 
                        tempdf=target[ selection ]
                        
                        try:
                            model=load(folder+'/'+ds+"/"+m+'/models/consolidated_'+df.iloc[i,1]+'_'+df.iloc[i,2]+'_model_trained.joblib')
                            preds = model.predict(tempdf) # predicting
                            tempdf['preds_binary']=preds
                                    
                            try:
                                preds_ = model.predict_proba(tempdf) # predicting probabilities
                                tempdf['preds_probs']=preds_
                            except:
                                pass  
                            name=df.iloc[i,1]+'_'+df.iloc[i,2]
                            tempdf.to_csv(folder+"/"+ds+"/featSelect_results/"+d+"/"+name+'_dataset_result.tsv',sep='\t')
                        except: 
                            pass
                
    def compile_results(self, folder_in, folder, controls):
        prefix=folder.split('/')[-1].split('_')[0]
        
        #folder="mtb_epitopes"    
        datasets_model={ 'epitope': [ 'bcipep','hla'] , 'protein': ['gram+', 'gram-','all_gram'] }
        methods=['method1','method2']
        passed=set()
        aux={}
        for d in datasets_model.keys():
            type_= d
            
            aux[d]={}
            f=open(f'{folder_in}/results/selected_{d}s.fasta','r')
            for line in f:
                l=line.replace('\n','')
                if(l.startswith('>')):
                    ide = l.replace('>','')
                    aux[d][ide]={}
            f.close()
            
            if(d=='epitope' and os.path.isfile(controls['control_epitope'])):
                aux[d]={}
                i=0
                f=open(controls['control_epitope'],'r')
                for line in f:
                    if(i>0):
                        l=line.replace('\n','').split('\t')
                        aux[d][l[0]]={'vaxijen_server': float(l[3]) }
                    i+=1
                f.close()
                
            if(d=='protein' and os.path.isfile(controls['control_protein']) ):
                aux[d]={}
                i=0
                f=open(controls['control_protein'],'r')
                for line in f:
                    if(i>0):
                        l=line.replace('\n','').split('\t')
                        aux[d][l[0]]={'vaxign-ml': float(l[1]) }
                    i+=1
                f.close()    
                
            for ds in datasets_model[d]:
                d=prefix+'_'+d.replace(prefix+'_','')
                
                for m in methods:
                    for res in os.listdir(folder+"/"+ds+"/best_models/results/"+d):
                        name = ds+'|'+m+'-'+res.replace('_dataset_result.tsv','')
                        df=pd.read_csv(folder+"/"+ds+"/best_models/results/"+d+"/"+res, sep='\t')
                        preds=df['preds_binary']
                        i=0
                        for k in aux[type_].keys():
                            aux[type_][k][name]=preds[i]
                            i+=1
                        
                    for res in os.listdir(folder+"/"+ds+"/featSelect_results/"+d):
                        name = ds+'|'+m+'-'+res.replace('_dataset_result.tsv','')
                        df=pd.read_csv(folder+"/"+ds+"/featSelect_results/"+d+"/"+res, sep='\t')
                        preds=df['preds_binary']
                        i=0
                        for k in aux[type_].keys():
                            aux[type_][k][name]=preds[i]
                            i+=1
        
        f=open(folder+'/compiled_predictions_antigenicity.tsv','w')
        f.write('sequence_type\tidentifier\tsource_dataset\tpredictor\tclass\n')
        for typ in aux.keys():
            for seqid in aux[typ].keys():
                for predictor in aux[typ][seqid].keys():
                    if(len(predictor.split('|'))==1):
                        source='external'
                        model=predictor
                    else:
                        source=predictor.split('|')[0]
                        model=predictor.split('|')[1]
                    f.write('%s\t%s\t%s\t%s\t%.6f\n' %(typ, seqid, source, model, aux[typ][seqid][predictor]) )
        f.close()
                
    def analyze_compiled(self, folder):
        prefix=folder.split('/')[-1].split('_')[0]
        #folder="mtb_epitopes"   
        flag = False 
        counts={}
        df=pd.read_csv(folder+'/compiled_predictions_antigenicity.tsv', sep='\t')
        for typ in df['sequence_type'].unique():
            counts[typ]={}      
            filt=df[ df['sequence_type']==typ ]
            total=len(filt['identifier'].unique())
            for ide in filt['identifier'].unique():
                ref=filt[ (filt['identifier']==ide) & (filt['source_dataset']=='external') ]
                if( len(ref) > 0 ):
                    ref = ref.iloc[0,-1]
                    dff=filt[ (filt['identifier']==ide) & (filt['source_dataset']!='external') ]
                    
                    for i in range(len(dff)):
                        source=dff.iloc[i, 2]
                        model=dff.iloc[i, 3]
                        
                        if(not source in counts[typ]):
                            counts[typ][source]={}
                        if(not model in counts[typ][source]):
                            counts[typ][source][model]=[0, total]
                        if(dff.iloc[i,-1]==ref):
                            counts[typ][source][model][0]+=1
                            
                        flag=True
        
        if(flag):             
            f=open(folder+'/checksum_predictions_antigenicity.tsv','w')
            f.write('sequence_type\tsource_dataset\tpredictor\tcount_agreement_ref\tpercentage_agreement\n')
            for typ in counts.keys():
                for source in counts[typ].keys():
                    for predictor in counts[typ][source].keys():
                        f.write('%s\t%s\t%s\t%i\t%.6f\n' %(typ, source, predictor, counts[typ][source][predictor][0], counts[typ][source][predictor][0]/counts[typ][source][predictor][1] ) )
            f.close()
    
    def get_statistics_external(self, ds):
        d=pd.read_csv(ds+'/compiled_predictions_antigenicity.tsv', sep='\t')
        pos=len( d[ (d['sequence_type']=='epitope') & (d['source_dataset']=='external') & (d['class']==1) ] )
        print('epitopes predicted positive: ', pos)
        pos=len( d[ (d['sequence_type']=='protein') & (d['source_dataset']=='external') & (d['class']==1) ] )
        print('proteins predicted positive: ', pos)
        
        proteins=d[ (d['sequence_type']=='protein')]['identifier'].unique()
        print('number of proteins: ', len(proteins) )
        epitopes=d[ (d['sequence_type']=='epitope')]['identifier'].unique()
        print('number of epitopes: ', len(epitopes) )
    
    def config_path_models(self, folder_in):
        dir_out=f"{folder_in}results/"
        
        if(folder_in.endswith('/')):
            folder_in = folder_in[:-1]
            
        prefix = folder_in.split('/')[-1].split('_')[0]
        if( not os.path.isdir( f"{folder_in}/{prefix}_antigenicity_prediction/") ):
            os.system( f"mkdir {folder_in}/{prefix}_antigenicity_prediction/" )
            
            os.system( f"mkdir {folder_in}/{prefix}_antigenicity_prediction/{prefix}_epitope" )
            os.system( f"cp {dir_out}selected_epitopes.fasta {folder_in}/{prefix}_antigenicity_prediction/{prefix}_epitope/dataset_epitope_test.fasta" )
            
            os.system( f"mkdir {folder_in}/{prefix}_antigenicity_prediction/{prefix}_protein" )
            os.system( f"cp {dir_out}selected_proteins.fasta {folder_in}/{prefix}_antigenicity_prediction/{prefix}_protein/dataset_protein_test.fasta" )
        
        os.system( f"cp -R trained_models/* {folder_in}/{prefix}_antigenicity_prediction/" )
    
    def get_consensus_antigenic_epitopes(self, ds, folder_in):
        dir_out=f"{folder_in}results/"
        
        cnt={}
        aux = pd.read_csv(ds+'/compiled_predictions_antigenicity.tsv', sep='\t')
        for i in aux.index:
            ide = aux.loc[i, 'identifier'] 
            c = aux.loc[i, 'class']  
            pr = aux.loc[i, 'predictor']  
            if(not ide in cnt):
                cnt[ide] = { 'c1': 0, 'cref': 0, 'ref': c }
                ant=ide
            else:
                if(c==1):
                    cnt[ide]['c1']+=1
                    
                    if(c==cnt[ide]['ref']):
                        cnt[ide]['cref']+=1
            
        external_protein = []
        external_epitope = []
        paprec_protein = []
        paprec_epitope = []
        agg_protein = []
        agg_epitope = []
        prots=set()
        res = pd.read_csv( f"{dir_out}table_final_proteins_epitopes.tsv", sep='\t')
        for i in res.index:
            ide = res.loc[i, 'epitope_id'] 
            idp = res.loc[i, 'protein'].split(',')[0]
            
            pape = cnt[ide]['c1'] > 1
            agge = cnt[ide]['cref'] > 1
            papp = cnt[idp]['c1'] > 1
            aggp = cnt[idp]['cref'] > 1
            
            """
            refe = aux[ (aux['identifier']==ide) & (aux['predictor']=='vaxijen_server') ].reset_index().loc[0, 'class']
            pape = len( aux[ (aux['identifier']==ide) & (aux['predictor']!='vaxijen_server') & (aux['class']==1) ] ) > 1
            agge = len( aux[ (aux['identifier']==ide) & (aux['predictor']!='vaxijen_server') & (aux['class']==refe) ] ) > 1
            
            
            refp = aux[ (aux['identifier']==idp) & (aux['predictor']=='vaxign-ml') ].reset_index().loc[0, 'class']
            papp = len( aux[ (aux['identifier']==idp) & (aux['predictor']!='vaxign-ml') & (aux['class']==1) ] ) > 1
            aggp = len( aux[ (aux['identifier']==idp) & (aux['predictor']!='vaxign-ml') & (aux['class']==refp) ] ) > 1
            """
            
            if( pape ):
                pape=1
            else:
                pape=0
                
            if( agge ):
                agge=1
            else:
                agge=0
                
            if( papp ):
                papp=1
            else:
                papp=0
                
            if( aggp ):
                aggp=1
            else:
                aggp=0
                
            paprec_epitope.append(pape)
            paprec_protein.append(papp)
            
            agg_epitope.append(agge)
            agg_protein.append(aggp)
            if(agge==1 and aggp==1):
                prots.add(idp)
                  
         
        res['agreement_protein_antigenic'] = agg_protein
        res['agreement_epitope_antigenic'] = agg_epitope
        res['paprec_protein_antigenic'] = paprec_protein
        res['paprec_epitope_antigenic'] = paprec_epitope
        cols = list(res.columns)[:-2]
        res = res[ (res['paprec_protein_antigenic'] == 1) | (res['paprec_epitope_antigenic'] == 1) ][ cols ]
        res.to_csv(f"{dir_out}table_final_antigenic_proteins_epitopes.tsv", sep='\t', index=None)
                
    def run(self):
        datasets = { "mtb_epitopes": { "control_protein": "epitocore/results_analysis/passed_proteins.result.tsv", "control_epitope": "epitocore/results_analysis/vaxijen2_prediction.tsv" }, "saureus_epitope":  { "control_protein": "epitocore/results_analysis_saureus/passed_proteins.result.tsv", "control_epitope": "epitocore/results_analysis_saureus/vaxijen2_prediction.tsv" } }
    
        if( self.config!=None ):
            for c in self.config['experiments']:
            
                folder_in = c['folder_in']
                folder = folder_in
                if(folder_in.endswith('/')):
                    folder = folder_in[:-1]
                else:
                    folder_in += '/'
                    
                c['control_protein'] = folder_in+'results/'+c['control_protein']
                c['control_epitope'] = folder_in+'results/'+c['control_epitope']
                
                prefix = folder.split('/')[-1].split('_')[0]
                d = f"{folder}/{prefix}_antigenicity_prediction"
                #if( not os.path.isdir(d) ):
                self.config_path_models( folder_in )
                
                print("---> Building datasets")
                self.build_exps_datasets(d) 
                
                print("---> Selecting models")
                self.select_models(d)
                
                print("---> Preparing test set")
                self.prepare_testeset(d)
                
                print("---> Using models to predict on new data separated")
                self.test_models_separated(d)
                
                print("---> Using models to predict on new data by feature selection")
                self.test_models_feature_selection(d)
                
                print("---> Compiling results")
                self.compile_results(folder, d, c)
                if( (c['control_protein'] != None and c['control_protein'] != '') or (c['control_epitope'] != None and c['control_epitope'] != '') ):
                    if( os.path.isfile(c['control_protein']) or os.path.isfile(c['control_epitope']) ):
                        #print("---> Analyzing compiled results")
                        #self.analyze_compiled(d)
                        
                        print("---> Comparing with the given reference prediction")
                        self.get_consensus_antigenic_epitopes(d, folder_in)
        else:
            print('Error: This module was not properly initialized')
        
    
import sys

if( len(sys.argv)==2 ):
    path_config = sys.argv[1]
    # step = int(sys.argv[2]) # 1 => filtering overlapping, 2 => parse il results
    a=AntigenicityTest( path_config)
    a.run()
else: 
    print('Error: incorrect number of passed parameters')
    
