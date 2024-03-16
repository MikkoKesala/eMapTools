from qgis.core import QgsVectorLayer,QgsField,edit
from qgis.PyQt.QtCore import QVariant
from .paras_2 import decay_tree_potential,NP_retention
from .geotools import hsAnalysis,go_gaussian,biodiversity
import numpy as np

def limitValue(x,minimum,maximum):
    """
    This adjust given value to minimum or maximum if outsite of the scale
    """
    if type(x)==float or type(x) == int:
        fc = lambda x : x if x>minimum and x<maximum else (minimum if x<=minimum else maximum)
    else:
        fc = lambda x : x
    l = fc(x)
    return l

def normalizeValue(in_feat:QgsVectorLayer,fieldname:str,filtervalues:tuple,transpose:bool=False):
    """
    This normalize values between 0 and 1 by formula f(x) = x-min(x) / max(x) - min(x)
    You can limit values between specific values and transpose the values
    """
    in_feat.dataProvider().addAttributes([QgsField(fieldname+"n",QVariant.Double)])
    in_feat.updateFields()
    
    lis = [feat[fieldname] for feat in in_feat.getFeatures() if feat[fieldname] is not None]
    
    #print (lis)
    if filtervalues is not None:
        lis = [limitValue(i,filtervalues[0],filtervalues[1]) for i in lis]
        with edit(in_feat):
            for feat in in_feat.getFeatures():
                feat[fieldname]=limitValue(feat[fieldname],filtervalues[0],filtervalues[1])
                in_feat.updateFeature(feat)
                
    try:
        minlis = min(list(lis))
        maxlis = max(list(lis))
        
    except:
        minlis = 0
        maxlis = 1
    
    
    with edit(in_feat):
        for feat in in_feat.getFeatures():
            if type(feat[fieldname]) in (int,float) and maxlis != minlis:
                if transpose == False:
                    feat[fieldname+"n"] = (feat[fieldname]-minlis) / (maxlis-minlis)
                else:
                    feat[fieldname+"n"] = 1- ((feat[fieldname]-minlis) / (maxlis-minlis))
            
            elif type(feat[fieldname]) in (int,float) and maxlis == minlis:
                feat[fieldname+"n"] = 0.0
            
    
            in_feat.updateFeature(feat)

def zScore(in_feat:QgsVectorLayer,fieldlist):
    """ This calculates z-scores every field in list """
    
    for field in fieldlist:
        in_feat.addAttributes([QgsField(field+"_z",QVariant.Double)])
        in_feat.updateFields()
    
    for field in fieldlist:
        vlist = [v[field] for v in in_feat.getFeatures() if v[field] is not None]
        s = np.std(vlist)
        m = np.mean(vlist)
        with edit(in_feat):
            for feat in in_feat.getFeatures():
                if feat[field] is not None:
                    feat[field+"_z"]=(m-feat[field]) / s
                
                in_feat.updateFeature(feat)

def treespeciesFromGrid(in_feat:QgsVectorLayer):
    """
    This determines treespecies from grid data. Grid data much have fields MEANDIAMETER<treespecies>
    """

    trees = {'MEANDIAMETERDECIDUOUS':3, 'MEANDIAMETERPINE':1, 'MEANDIAMETERSPRUCE':2}

    
    in_feat.dataProvider().addAttributes([QgsField("treespecies",QVariant.Int)])
    in_feat.dataProvider().addAttributes([QgsField("diameter",QVariant.Double)])
    in_feat.updateFields()
    
    with edit(in_feat):
        for c,feat in enumerate(in_feat.getFeatures()):
            m = max([feat[t] for t in trees])
            #print (m)
            i  = [t for t in trees if feat[t]==m]
            #print (m,i,trees[i[0]])
            feat["treespecies"]=trees[i[0]]
            feat["diameter"]=m
            #feat['biod'] = float(sim_di+conver_cof)
            in_feat.updateFeature(feat)

def treespeciesFromGrid2(in_feat:QgsVectorLayer,chm_height):
    trees = {'MEANDIAMETERDECIDUOUS':3, 'MEANDIAMETERPINE':1, 'MEANDIAMETERSPRUCE':2}

    u = 10

    in_feat.dataProvider().addAttributes([QgsField("treespecies",QVariant.Int)])
    in_feat.dataProvider().addAttributes([QgsField("diameter",QVariant.Double)])
    in_feat.dataProvider().addAttributes([QgsField("chm_height",QVariant.Double)])
    in_feat.updateFields()
    
    with edit(in_feat):
        for c,feat in enumerate(in_feat.getFeatures()):
            try:
                m = max([feat[t] for t in trees if type(feat[t]) in (int,float)])
            except:
                m = None
            h_list = [round(abs(feat[chm_height]/u-feat[t.replace('DIAMETER','HEIGHT')]),2) for t in trees if type(feat[t.replace('DIAMETER','HEIGHT')]) in (int,float) and feat[t]>0] #getting closest difference value between chm_height and grid heights
            if len(h_list)>0:
                h = min(h_list)
                i =  [t for t in trees if type(feat[t.replace('DIAMETER','HEIGHT')]) in (int,float) and round(abs(feat[chm_height]/u-feat[t.replace('DIAMETER','HEIGHT')]),2) == h] #getting dict key of the closest height value 
                feat["treespecies"]=trees[i[0]]
                feat["diameter"]=(feat[chm_height]/u) / feat[i[0].replace('DIAMETER','HEIGHT')]*feat[i[0]] #diameter = chm_height / meanheight<treepspecies> * meandiameter<treepspecies>
                feat['chm_height'] = feat[chm_height]/u
            elif type(m) in (float,int):
                #m = max([feat[t] for t in trees])
                i  = [t for t in trees if feat[t]==m]
                feat["treespecies"]=trees[i[0]]
                feat["diameter"]=m
                feat['chm_height'] = feat[chm_height]/u

            
            in_feat.updateFeature(feat)

def simpson_di(species):
    """
    This calculates simpson diversity index by formula D = n(n-1)/N(N-1)
    Species list are form of [10,1,3,4]
    """
    
    proportions = [i/sum(species) for i in species]
    simpsons_index = 0
    for proportion in proportions:
        simpsons_index += proportion ** 2
    
    return (1 - simpsons_index) / (1-1/len(species)) if simpsons_index > 0 and len(species)>1 else 0

def calculateBiodiversity(in_feat:QgsVectorLayer,speciesfield:list):
    """
    This calculates specific diversity index by given species list and save it to given QgsVectorLayer.
    Species list have to be format of [1,35,42,2,23,...,n]
    """
    
    in_feat.dataProvider().addAttributes([QgsField("biod",QVariant.Double),QgsField("treedistribution",QVariant.List)])
    in_feat.updateFields()
    
    with edit(in_feat):
        for feat in in_feat.getFeatures():
            
            sim_di = simpson_di([feat[i] for i in speciesfield if type(feat[i]) in (int,float)])
            dis = [feat[i] for i in speciesfield if type(feat[i]) in (int,float)]
            dis = [round(i / sum(dis),2) for i in dis]
            feat['biod'] = round(float(sim_di),3)
            feat['treedistribution']=dis
            in_feat.updateFeature(feat)
    
    #normalizeValue(in_feat,'biod',None,False)

def calculateDecayTreePotential(in_feat,fz_field):
    
    in_feat.dataProvider().addAttributes([QgsField("dtree",QVariant.Double)])
    in_feat.updateFields()
    
    #treelist = [1,2,29]
    puuH = [("MEANHEIGHTPINE",1),("MEANHEIGHTSPRUCE",2),("MEANHEIGHTDECIDUOUS",29)]
    fc_d = lambda t : t.replace("HEIGHT","DIAMETER")
    fc_v = lambda t : t.replace("HEIGHT","VOLUME")
    
    with edit(in_feat):
        for feat in in_feat.getFeatures():
            if type(feat[fz_field]) in (str,int,float):
                dcp = decay_tree_potential('zone'+str(feat[fz_field]))
            else:
                dcp = decay_tree_potential('zone3')
            if type(feat['FERTILITYCLASS']) in (float,int) and max([feat[p[0]] for p in puuH]) >0:
                if feat['FERTILITYCLASS']>6:
                    para = [dcp[6][i[1]] for i in puuH]
                else:
                    para = [dcp[int(feat['FERTILITYCLASS'])][i[1]] for i in puuH]
                                    
                potvalues = [limitValue(np.poly1d(para[i])(feat[fc_d(p[0])]),0,2) for i,p in enumerate(puuH) if feat[fc_d(p[0])] > 0]
            
            else:
                potvalues=[0]
            
            feat["dtree"]=float(sum(potvalues))
            
            
            in_feat.updateFeature(feat)
    
    #normalizeValue(in_feat,"dtree",None,False)

def decay2tree(in_feat:QgsVectorLayer,diameter:str,fertilityclass:str,treespecies:str,biogeoclass:str):
    """
    This calcultes decaytree value for single tree. Input layer is in qgsvectorlayer format and other parameteres are fieldnames.
    """
    in_feat.dataProvider().addAttributes([QgsField("dtree",QVariant.Double)])
    in_feat.updateFields()
    

    with edit(in_feat):
        for feat in in_feat.getFeatures():
            if type(feat[biogeoclass]) in (str,int,float):
                dcp = decay_tree_potential('zone'+str(feat[biogeoclass]))
            else:
                dcp = decay_tree_potential('zone3')
            #dcp = decay_tree_potential('zone'+str(feat[biogeoclass]))

            if type(feat[fertilityclass]) in (int,float) and feat[diameter] > 0 and type(feat[treespecies]) in (int,float):
                if feat[fertilityclass]>6:
                    para = dcp[6][feat[treespecies]]
                else:
                    para = dcp[feat[fertilityclass]][feat[treespecies]]
                                    
                potvalues = np.poly1d(para)(feat[diameter])
                #potvalues = limitValue(potvalues,0,2)
            else:
                potvalues=0
            
            feat["dtree"]=round(float(potvalues),3)
            
            
            in_feat.updateFeature(feat)
    
    #normalizeValue(in_feat,"dtree",None,False)


def calculateNPretention(in_feat):
    in_feat.dataProvider().addAttributes([QgsField("pRetent",QVariant.Double)])
    in_feat.updateFields()
    ret = NP_retention()
    min_dist = np.min([i['wdistance'] for i in in_feat.getFeatures()])
    if min_dist <= 50:
        with edit(in_feat):
            for feat in in_feat.getFeatures():
                if type(feat['wdistance']) in (float,int) and feat['wdistance']>1:
                    feat['pRetent'] = ret['P']/feat['wdistance']
                elif feat['wdistance']==0:
                    feat['pRetent']=ret['P']
                
                in_feat.updateFeature(feat)
        
    #normalizeValue(in_feat,'pRetent',None,False)

def calculateEnvValue(in_feat,weights):
    in_feat.dataProvider().addAttributes([QgsField("env_value",QVariant.Double)])
    in_feat.updateFields()
    ecovalues = ['biodn','pRetentn','DTW_1n','dtreen']
    with edit(in_feat):
        for feat in in_feat.getFeatures():
            eco = []
            for e in ecovalues:
                if type(feat[e]) in (int,float):
                    eco.append(feat[e])
                else:
                    eco.append(0)
            feat['env_value'] = eco[0]*weights['BIO']+eco[1]*weights['NP']+eco[2]*weights['DTW']+eco[3]*weights['LP']
            
            #normalizeValue(in_feat,'env_value',None,False)
            in_feat.updateFeature(feat)

def selectReTrees(in_feat:QgsVectorLayer,fieldname:str,cuttingfield:str,treecount:int,cuttingsize:float):
    """
    This select retention trees by given QgsVectorLayer, fieldname and treecount. Parameter "cuttingfield" restrict selection to cutting area 
    """
    
    in_feat.dataProvider().addAttributes([QgsField("reTree",QVariant.Int)])
    in_feat.updateFields()

    #Restrict to cutting area
    NoLeim = [feat.id() for feat in in_feat.getFeatures() if feat[cuttingfield]!=1]
    in_feat.dataProvider().deleteFeatures(NoLeim)
    in_feat.updateFields()

    #getting to values in which we do the selection
    opt = np.array([feat[fieldname] for feat in in_feat.getFeatures() if type(feat[fieldname]) in (int,float)])

    treecount = int(round(float(treecount) * float(cuttingsize),0))
    pvalue = opt[np.argsort(opt)[-treecount]]
    
    
    with edit(in_feat):
        for feat in in_feat.getFeatures():
            if feat[fieldname]>=pvalue:
                feat['reTree']=1
            else:
                feat['reTree']=0
            in_feat.updateFeature(feat)

def runEssModel(in_feat:QgsVectorLayer,weights,treecount,cuttingsize,fz_field):
    normalizeValue(in_feat,"DTW_1",(0.0,1.0),True)
    treespeciesFromGrid2(in_feat,"CHM")
    in_feat = biodiversity(in_feat,'fid','treespecies',20)
    #calculateBiodiversity(in_feat,["STEMCOUNTPINE","STEMCOUNTDECIDUOUS","STEMCOUNTSPRUCE"])
    decay2tree(in_feat,'diameter','fertilityclass','treespecies',fz_field)
    #calculateDecayTreePotential(in_feat,fz_field)
    calculateNPretention(in_feat)
    normalizeValue(in_feat,"biod",None,False)
    normalizeValue(in_feat,"dtree",(0.0,2.0),False)
    normalizeValue(in_feat,"pRetent",None,False)
    calculateEnvValue(in_feat,weights)
    #retrees = hsAnalysis(in_feat,'env_value')
    #id_name = in_feat.fields().id()
    retrees = go_gaussian(in_feat,'fid','env_value',30)
    #normalizeValue(retrees,"getisord",None,False)
    selectReTrees(retrees,'getisord','leimikko',treecount,cuttingsize)

    return retrees

def runEssModel2points(in_feat:QgsVectorLayer,weights,treecount,cuttingsize,fz_field,d_field,s_field):
    normalizeValue(in_feat,"DTW_1",(0.0,1.0),True)
    #treespeciesFromGrid2(in_feat,"CHM")
    in_feat = biodiversity(in_feat,'fid',s_field,20)
    #calculateBiodiversity(in_feat,["STEMCOUNTPINE","STEMCOUNTDECIDUOUS","STEMCOUNTSPRUCE"])
    decay2tree(in_feat,d_field,'fertilityclass',s_field,fz_field)
    #calculateDecayTreePotential(in_feat,fz_field)
    calculateNPretention(in_feat)
    normalizeValue(in_feat,"biod",None,False)
    normalizeValue(in_feat,"dtree",(0.0,2.0),False)
    normalizeValue(in_feat,"pRetent",None,False)
    calculateEnvValue(in_feat,weights)
    #retrees = hsAnalysis(in_feat,'env_value')
    #id_name = in_feat.fields().id()
    retrees = go_gaussian(in_feat,'fid','env_value',30)
    #normalizeValue(retrees,"getisord",None,False)
    selectReTrees(retrees,'getisord','leimikko',treecount,cuttingsize)
    return retrees

def runEssModel2points2(in_feat:QgsVectorLayer,weights,cuttingsize,attributes):
    #attnames = ['species','diameter','treecount','stemcounts','vegetationzones','fertilityclass']
    normalizeValue(in_feat,"DTW_1",(0.0,1.0),True)
    #treespeciesFromGrid2(in_feat,"CHM")
    calculateBiodiversity(in_feat,attributes['stemcounts'])
    decay2tree(in_feat,attributes['diameter'],attributes['fertilityclass'],attributes['species'],attributes['vegetationzones'])
    calculateNPretention(in_feat)
    normalizeValue(in_feat,"biod",None,False)
    normalizeValue(in_feat,"dtree",None,False)
    normalizeValue(in_feat,"pRetent",None,False)
    calculateEnvValue(in_feat,weights)
    retrees = hsAnalysis(in_feat,'env_value')
    normalizeValue(retrees,"HS_1",None,False)
    selectReTrees(retrees,'HS_1n','leimikko',attributes['treecount'],cuttingsize)

    return retrees