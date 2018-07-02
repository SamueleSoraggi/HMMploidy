
#! /usr/bin/python3

''' file to infer the critical values of delta likelihood to infer
 aneuploidy '''

import pandas as pd
import re
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
from sklearn.preprocessing import OneHotEncoder

#import the data
Data=pd.read_csv('Interpret_delta_new',sep='\t',names=["Ploidies","Inferred_Ploidy","SNPs","Mean_read_depth","H_0","H_1","Delta","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"],index_col=False)
# Create new column with normalised data and other features wanted
Data['Normalised_delta']=(Data['Delta']/Data['SNPs'])
Ploids=list(Data['Ploidies'])
base=[]
no_base=[]
aneu=[]
no_aneu=[]
is_aneu=[]
NSAMS=[]
for ploid in Ploids:
    ploid=re.split(',|x',ploid)
    base.append(int(ploid[0]))
    no_base.append(int(ploid[1]))
    aneu.append(int(ploid[2]))
    no_aneu.append(int(ploid[3]))
    if int(ploid[0])==int(ploid[2]):
        is_aneu.append(0)
    else:
        is_aneu.append(1)
    NSAMS.append(int(ploid[1])+int(ploid[3]))

Data['Base_ploidy']=base
Data['Base_ploidy_number']=no_base
Data['Aneuploidy_level']=aneu
Data['Number_of_aneuploidy']=no_aneu
Data['NSAMS']=NSAMS
Data['Aneuploidy']=is_aneu




#Seperate Data into X: the independent variables and y the dependent variables
X_All=Data[['Inferred_Ploidy','Mean_read_depth','Normalised_delta','NSAMS','H_0','H_1','SNPs']].values
X=Data[['Inferred_Ploidy','Mean_read_depth','Normalised_delta','NSAMS']].values
y=Data.iloc[:,-1].values
#X=X[1800:]
#y=y[1800:]

# Encode categorical data
#from sklearn.preprocessing import OneHotEncoder
#onehotencoder=OneHotEncoder(categorical_features=[0])
#X=onehotencoder.fit_transform(X).toarray()

from sklearn.preprocessing import PolynomialFeatures 
polynomial_regressor = PolynomialFeatures(degree=2) # specify want terms up to degree 2
X_poly = polynomial_regressor.fit_transform(X)  # transform data

# Seperate data into a training set and set to test on
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X_poly,y, test_size=0.25, random_state=1)


# Fit a logistic regression model to the training set
from sklearn.linear_model import LogisticRegression
classifier = LogisticRegression(random_state=0)
classifier.fit(X_train,y_train)

#Predict the test set results
y_pred=classifier.predict(X_test)

#Evaluating the performance of the logistic regression using confusion matrix
from sklearn.metrics import confusion_matrix
cm=confusion_matrix(y_test,y_pred)
print(cm)
accuracy=(cm[0][0]+cm[1][1])/sum(sum(cm))
print(accuracy)  
#can get predictiosn using classifier.predict([[data]])

import pickle
pickle.dump(classifier, open('Aneuploidy_Classifier.sav', 'wb'))
pickle.dump(polynomial_regressor, open('Polynomial_Regressor.pk1','wb'))
 
