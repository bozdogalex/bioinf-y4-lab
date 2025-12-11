**Demo01:**
[Supervised] Test accuracy: 0.944

Crosstab TrueLabel vs Cluster:
Cluster     0   1   2
TrueLabel            
0           0  98   2
1           2  11  87
2          82  18   0


**Demo02:**
[INFO] Loading demo data from labs/08_ML_flower/submissions/codmihaic/data/expression_demo.csv
      Gene_1     Gene_2     Gene_3     Gene_4     Gene_5     Gene_6     Gene_7  ...    Gene_75    Gene_76    Gene_77    Gene_78    Gene_79    Gene_80  Label
0  12.496714  11.861736  12.647689  13.523030  11.765847  11.765863  13.579213  ...   5.380255   8.821903   8.087047   7.700993   8.091761   6.012431  Brain
1   7.780328   8.357113   9.477894   7.481730   7.191506   7.498243   8.915402  ...   8.293072   7.285649   9.865775   8.473833   6.808697   8.656554  Liver
2   7.025318   8.787085   9.158596   7.179318   8.963376   8.412781   8.822060  ...  10.143944   8.633919   5.974857   8.186454   7.338214   8.852433  Heart
3   7.207479   7.885264   8.504987   8.865755   6.799704   7.665499   7.525055  ...  13.305479  12.021004  12.681953  11.689733  12.324166  11.869857  Tumor
4  12.096996  12.595157  11.181779  14.092387  10.993983  10.785811  13.158111  ...   9.179440   7.530824   6.286865   9.353872   7.885460   9.237816  Brain

[5 rows x 81 columns]
[INFO] Classes: ['Brain', 'Heart', 'Liver', 'Tumor']

=== Classification Report (RF) ===
              precision    recall  f1-score   support

       Brain       1.00      1.00      1.00        12
       Heart       1.00      1.00      1.00        12
       Liver       1.00      1.00      1.00        12
       Tumor       1.00      1.00      1.00        12

    accuracy                           1.00        48
   macro avg       1.00      1.00      1.00        48
weighted avg       1.00      1.00      1.00        48


Crosstab Label vs Cluster:
Cluster   0   1   2   3
Label                  
Brain     0  60   0   0
Heart     0   0  60   0
Liver    60   0   0   0
Tumor     0   0   0  60
[INFO] Demo finished.


**Demo03:**
[INFO] Loading demo data from labs/08_ML_flower/submissions/codmihaic/data/expression_demo.csv
      Gene_1     Gene_2     Gene_3     Gene_4     Gene_5     Gene_6     Gene_7  ...    Gene_75    Gene_76    Gene_77    Gene_78    Gene_79    Gene_80  Label
0  12.496714  11.861736  12.647689  13.523030  11.765847  11.765863  13.579213  ...   5.380255   8.821903   8.087047   7.700993   8.091761   6.012431  Brain
1   7.780328   8.357113   9.477894   7.481730   7.191506   7.498243   8.915402  ...   8.293072   7.285649   9.865775   8.473833   6.808697   8.656554  Liver
2   7.025318   8.787085   9.158596   7.179318   8.963376   8.412781   8.822060  ...  10.143944   8.633919   5.974857   8.186454   7.338214   8.852433  Heart
3   7.207479   7.885264   8.504987   8.865755   6.799704   7.665499   7.525055  ...  13.305479  12.021004  12.681953  11.689733  12.324166  11.869857  Tumor
4  12.096996  12.595157  11.181779  14.092387  10.993983  10.785811  13.158111  ...   9.179440   7.530824   6.286865   9.353872   7.885460   9.237816  Brain

[5 rows x 81 columns]
[INFO] Classes: ['Brain', 'Heart', 'Liver', 'Tumor']
/usr/local/lib/python3.11/site-packages/sklearn/linear_model/_logistic.py:1247: FutureWarning: 'multi_class' was deprecated in version 1.5 and will be removed in 1.7. From then on, it will always use 'multinomial'. Leave it to its default value to avoid this warning.
  warnings.warn(

=== Classification Report (Logistic Regression) ===
              precision    recall  f1-score   support

       Brain       1.00      1.00      1.00        12
       Heart       1.00      1.00      1.00        12
       Liver       1.00      1.00      1.00        12
       Tumor       1.00      1.00      1.00        12

    accuracy                           1.00        48
   macro avg       1.00      1.00      1.00        48
weighted avg       1.00      1.00      1.00        48

[INFO] Saved confusion matrix to labs/08_ML_flower/submissions/codmihaic/demo_outputs/demo_logreg_confusion.png
[INFO] Logistic Regression demo finished.


**Task01:**
              precision    recall  f1-score   support

       Brain       0.88      0.96      0.92        23
       Heart       0.67      0.29      0.40         7
       Liver       0.70      0.82      0.76        17
       Tumor       1.00      0.80      0.89         5

    accuracy                           0.81        52
   macro avg       0.81      0.72      0.74        52
weighted avg       0.80      0.81      0.79        52

Crosstab Label vs Cluster:
Cluster   0   1   2   3
Label                  
Brain    39  66   1  10
Heart     1   0  23  12
Liver    10  25  14  37
Tumor     0   0  20   2


**Task02:**
/usr/local/lib/python3.11/site-packages/sklearn/linear_model/_logistic.py:1247: FutureWarning: 'multi_class' was deprecated in version 1.5 and will be removed in 1.7. From then on, it will always use 'multinomial'. Leave it to its default value to avoid this warning.
  warnings.warn(
=== Random Forest ===
              precision    recall  f1-score   support

       Brain       0.88      0.96      0.92        23
       Heart       0.67      0.29      0.40         7
       Liver       0.70      0.82      0.76        17
       Tumor       1.00      0.80      0.89         5

    accuracy                           0.81        52
   macro avg       0.81      0.72      0.74        52
weighted avg       0.80      0.81      0.79        52


=== Logistic Regression ===
              precision    recall  f1-score   support

       Brain       0.96      0.96      0.96        23
       Heart       0.43      0.43      0.43         7
       Liver       0.75      0.88      0.81        17
       Tumor       1.00      0.40      0.57         5

    accuracy                           0.81        52
   macro avg       0.78      0.67      0.69        52
weighted avg       0.82      0.81      0.80        52