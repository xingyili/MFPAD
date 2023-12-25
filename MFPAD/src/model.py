import pandas as pd
import xgboost as xgb
import numpy as np
from xgboost import XGBRegressor
from sklearn.metrics import mean_squared_error, make_scorer, confusion_matrix
from sklearn.model_selection import KFold
import math
x=pd.read_excel("Feature.xlsx")
y=pd.read_excel("distance.xlsx")
accuracy=[]

#1. 划分
kf = KFold(n_splits=10, random_state=2023, shuffle=True)
accuracy_0 = []
accuracy_1 = []
acc_all=[]
rmse=[]
sen_all=[]
spe_all=[]
acc_2022_2023_all=[]
feature_importance_list=[]
for index, (train_index, test_index) in enumerate(kf.split(x, y)):
    x_train, x_test = np.array(x)[train_index], np.array(x)[test_index]
    y_train, y_test = np.array(y)[train_index], np.array(y)[test_index]
    xgb_model = xgb.XGBRegressor(max_depth=13,
                             booster='gbtree',
                             gamma='0.1',
                            learning_rate=0.1,
                            n_estimators=200,
                            objective='reg:squarederror',
                             seed=250)


    xgb_model.fit(x_train,y_train)
    feature_importance_list.append(xgb_model.feature_importances_)
    y_test_pre = xgb_model.predict(x_test)
    msetest = mean_squared_error(y_test, y_test_pre)
    rmse.append(msetest)
    print(msetest)
    # array = pd.to_numeric(y_test["distance"])  # object类型转为int类型
    # y_test = array.tolist()
    y_test_b = []
    # print(len(y_test))
    for i in y_test:
        if i > 2:
            flag = 1
        else:
            flag = 0
        y_test_b.append(flag)
    y_pre_b = []
    # print(len(y_test_pre))
    for i in y_test_pre:
        if i> 2:
            flag = 1
        else:
            flag = 0
        y_pre_b.append(flag)
    acc_1 = 0
    acc_0 = 0
    len_0=0
    len_1=0
    for i in range(len(y_test_b)):
        if y_test_b[i]==0:
            len_0=len_0+1
        if y_test_b[i]==1:
            len_1=len_1+1
        if y_test_b[i] == y_pre_b[i] and y_test_b[i]==1:
            acc_1= acc_1 + 1
        if y_test_b[i] == y_pre_b[i] and y_test_b[i]==0:
            acc_0= acc_0 + 1
    accuracy_0.append(acc_0/len_0)
    accuracy_1.append(acc_1/len_1)
    print(f"acc_0:{acc_0/len_0}")
    # print(acc_1/len_1)
    print(f"acc_1:{acc_1/len_1}")
    print(f"acc_1:{acc_1+acc_0 / len(y_test_pre)}")
    # acc_all.append((acc_0+acc_1)/len(y_test_pre))

    conf_matrix = confusion_matrix(y_test_b, y_pre_b)

    sen = conf_matrix[1, 1] / (conf_matrix[1, 1] + conf_matrix[1, 0])
    spe = conf_matrix[0, 0] / (conf_matrix[0, 0] + conf_matrix[0, 1])
    sen_all.append(sen)
    spe_all.append(spe)
    # 打印结果
    print(f"Sensitivity: {sen}")
    print(f"Specificity: {spe}")

for i in range(len(rmse)):
    rmse[i]=math.sqrt(rmse[i])

print("all:")
print("accuracy_0:")
print(sum(accuracy_0)/10)
print("accuracy_1:")
print(sum(accuracy_1)/10)
print("acc_all:")
print(sum(acc_all)/10)
print("rmse_all:")
print(sum(rmse)/10)
print("sen_all:")
print(sum(sen_all)/10)
print("spe_all:")
print(sum(spe_all)/10)








