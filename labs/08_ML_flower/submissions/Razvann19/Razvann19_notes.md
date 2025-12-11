**DEMO 1**
[Supervised] Test accuracy: 0.944

Crosstab TrueLabel vs Cluster:
Cluster     0   1   2
TrueLabel            
0           0  98   2
1           2  11  87
2          82  18   0

**DEMO 2**
[INFO] Loading demo data from labs/08_ML_flower/submissions/Razvann19/expression_demo.csv
     G1    G2    G3    G4    G5    G6    G7    G8    G9   G10    Label
0  5.11  3.86  6.13  4.70  2.36  7.28  6.18  7.81  2.03  5.52  TissueA
1  5.10  3.99  6.04  4.33  2.67  7.02  6.32  8.11  2.37  5.33  TissueA
2  4.88  4.07  6.28  4.70  2.54  6.92  6.15  7.98  2.31  5.46  TissueA
3  4.75  4.23  5.98  4.44  2.15  6.75  6.42  8.01  2.11  5.11  TissueA
4  4.95  4.25  5.96  4.27  2.26  6.96  5.90  7.99  2.27  5.40  TissueA
[INFO] Classes: ['TissueA', 'TissueB', 'TissueC']

=== Classification Report (RF) ===
              precision    recall  f1-score   support

     TissueA       1.00      1.00      1.00         2
     TissueB       1.00      1.00      1.00         2
     TissueC       1.00      1.00      1.00         2

    accuracy                           1.00         6
   macro avg       1.00      1.00      1.00         6
weighted avg       1.00      1.00      1.00         6


Crosstab Label vs Cluster:
Cluster   0   1   2
Label              
TissueA   0   0  10
TissueB  10   0   0
TissueC   0  10   0
[INFO] Demo finished.


**DEMO 3**


[INFO] Loading demo data from labs/08_ML_flower/submissions/Razvann19/expression_demo.csv
     G1    G2    G3    G4    G5    G6    G7    G8    G9   G10    Label
0  5.11  3.86  6.13  4.70  2.36  7.28  6.18  7.81  2.03  5.52  TissueA
1  5.10  3.99  6.04  4.33  2.67  7.02  6.32  8.11  2.37  5.33  TissueA
2  4.88  4.07  6.28  4.70  2.54  6.92  6.15  7.98  2.31  5.46  TissueA
3  4.75  4.23  5.98  4.44  2.15  6.75  6.42  8.01  2.11  5.11  TissueA
4  4.95  4.25  5.96  4.27  2.26  6.96  5.90  7.99  2.27  5.40  TissueA
[INFO] Classes: ['TissueA', 'TissueB', 'TissueC']
/usr/local/lib/python3.11/site-packages/sklearn/linear_model/_logistic.py:1247: FutureWarning: 'multi_class' was deprecated in version 1.5 and will be removed in 1.7. From then on, it will always use 'multinomial'. Leave it to its default value to avoid this warning.
  warnings.warn(

=== Classification Report (Logistic Regression) ===
              precision    recall  f1-score   support

     TissueA       1.00      1.00      1.00         2
     TissueB       1.00      1.00      1.00         2
     TissueC       1.00      1.00      1.00         2

    accuracy                           1.00         6
   macro avg       1.00      1.00      1.00         6
weighted avg       1.00      1.00      1.00         6

[INFO] Saved confusion matrix to labs/08_ml/demo_outputs/demo_logreg_confusion.png
[INFO] Logistic Regression demo finished.


** EX 01**

[INFO] Găsit fișierul de input: labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv
[INFO] Încărcăm dataset-ul din labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv
[INFO] Shape dataframe: (1000, 22)
[INFO] Coloane: ['Unnamed: 0', 'Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5', 'Sample_6', 'Sample_7', 'Sample_8', 'Sample_9', 'Sample_10', 'Sample_11', 'Sample_12', 'Sample_13', 'Sample_14', 'Sample_15', 'Sample_16', 'Sample_17', 'Sample_18', 'Sample_19', 'Sample_20', 'Label']
[INFO] Număr features numerice: 20
[INFO] Primele etichete unice: ['medium' 'high' 'low']
[INFO] Clase encodate: ['high', 'low', 'medium']
[INFO] Antrenăm RandomForestClassifier...
[INFO] Antrenare RF finalizată.
[INFO] Evaluăm modelul pe setul de test...

=== Classification Report (Random Forest) ===
              precision    recall  f1-score   support

        high       1.00      1.00      1.00        67
         low       1.00      1.00      1.00        67
      medium       1.00      1.00      1.00        66

    accuracy                           1.00       200
   macro avg       1.00      1.00      1.00       200
weighted avg       1.00      1.00      1.00       200

[INFO] Salvat classification_report în labs/08_ML_flower/submissions/Razvann19/classification_report_Razvann19.txt
[INFO] Salvat matricea de confuzie în labs/08_ML_flower/submissions/Razvann19/confusion_rf_Razvann19.png
[INFO] Calculăm importanța trăsăturilor (feature_importances_)...
[INFO] Salvat feature importance în labs/08_ML_flower/submissions/Razvann19/feature_importance_Razvann19.csv

[INFO] Top 20 features după importanță:
      Feature  Importance
19  Sample_20    0.698314
6    Sample_7    0.018200
3    Sample_4    0.017404
9   Sample_10    0.017286
18  Sample_19    0.017056
10  Sample_11    0.017028
14  Sample_15    0.016771
13  Sample_14    0.016573
8    Sample_9    0.016263
11  Sample_12    0.016183
2    Sample_3    0.015936
1    Sample_2    0.015557
7    Sample_8    0.015534
4    Sample_5    0.015084
5    Sample_6    0.014916
15  Sample_16    0.014843
12  Sample_13    0.014601
0    Sample_1    0.014544
16  Sample_17    0.014067
17  Sample_18    0.013839
[INFO] Rulăm KMeans cu n_clusters=3...
[INFO] Crosstab Label vs Cluster:
Cluster    0  1    2
Label               
high     204  0  129
low      287  0   47
medium   273  1   59
[INFO] Salvat crosstab-ul în labs/08_ML_flower/submissions/Razvann19/cluster_crosstab_Razvann19.csv
[INFO] Exercise 8 terminat cu succes.

** EX 02 **

[INFO] Găsit fișierul de input: labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv
[INFO] Încărcăm dataset-ul din labs/08_ML_flower/submissions/Razvann19/expression_matrix.csv
[INFO] Shape dataframe: (1000, 22)
[INFO] Coloane: ['Unnamed: 0', 'Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5', 'Sample_6', 'Sample_7', 'Sample_8', 'Sample_9', 'Sample_10', 'Sample_11', 'Sample_12', 'Sample_13', 'Sample_14', 'Sample_15', 'Sample_16', 'Sample_17', 'Sample_18', 'Sample_19', 'Sample_20', 'Label']
[INFO] Număr features numerice: 20
[INFO] Primele etichete unice: ['medium' 'high' 'low']
[INFO] Clase encodate: ['high', 'low', 'medium']
[INFO] Antrenăm modelele RF și Logistic Regression...
/usr/local/lib/python3.11/site-packages/sklearn/linear_model/_logistic.py:1247: FutureWarning: 'multi_class' was deprecated in version 1.5 and will be removed in 1.7. From then on, it will always use 'multinomial'. Leave it to its default value to avoid this warning.
  warnings.warn(
[INFO] Antrenare modele finalizată.
[INFO] Comparam modelele pe setul de test...
=== Random Forest ===
              precision    recall  f1-score   support

        high       1.00      1.00      1.00        67
         low       1.00      1.00      1.00        67
      medium       1.00      1.00      1.00        66

    accuracy                           1.00       200
   macro avg       1.00      1.00      1.00       200
weighted avg       1.00      1.00      1.00       200


=== Logistic Regression ===
              precision    recall  f1-score   support

        high       0.98      0.96      0.97        67
         low       0.93      1.00      0.96        67
      medium       0.95      0.91      0.93        66

    accuracy                           0.95       200
   macro avg       0.96      0.95      0.95       200
weighted avg       0.96      0.95      0.95       200

[INFO] Salvat raport comparativ în labs/08_ML_flower/submissions/Razvann19/rf_vs_logreg_report_Razvann19.txt
[INFO] Exercise 8b — RF vs Logistic Regression rulat cu succes.