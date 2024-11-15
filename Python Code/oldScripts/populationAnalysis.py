
fileList = glob.glob('../Meetz/corrMaster/*.npy')
masterMat = np.load(fileList[0])

for i in fileList[1:]:
    tempMat = np.load(i)
    masterMat = np.concatenate((masterMat,tempMat),axis=1)


## combined average score
plt.scatter(combCombinedSimScore, combCorrMat, color='green')
# regression line
m, b = np.polyfit(combCombinedSimScore, combCorrMat, 1)
#add linear regression line to scatterplot 
plt.plot(combCombinedSimScore, m*combCombinedSimScore+b)
plt.xlabel('Population Average Sim Score (dir tuning + RF location)')
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.ylabel('Population Pearson Correlation')
plt.xlim([0,1])
pearsonR, pValue = stats.pearsonr(combCombinedSimScore,combCorrMat)
plt.axis('equal')
plt.title(f'Pearson Correalation of X/Y axis = {pearsonR:.2f}, p-value = {pValue:.2g}', fontsize=8)
plt.show()


## combined average score (pref dir and RF overlap)
plt.scatter(combCombinedSimScorePref, combCorrMat, color='green')
# regression line
m, b = np.polyfit(combCombinedSimScorePref, combCorrMat, 1)
#add linear regression line to scatterplot 
plt.plot(combCombinedSimScorePref, m*combCombinedSimScorePref+b)
plt.xlabel('Population Average Sim Score (Pref Dir + RF location)')
plt.ylabel('Population Pearson Correlation')
plt.xlim([0,1])
pearsonR, pValue = stats.pearsonr(combCombinedSimScorePref,combCorrMat)
plt.axis('equal')
plt.title(f'Pearson Correalation of X/Y axis = {pearsonR}, p-value = {pValue}', fontsize=8)
plt.show()


## combined Multiplied combined sim score 
plt.scatter(combCombinedSimScoreMulti, combCorrMat, color='green')
# regression line
m, b = np.polyfit(combCombinedSimScoreMulti, combCorrMat, 1)
#add linear regression line to scatterplot 
plt.plot(combCombinedSimScoreMulti, m*combCombinedSimScoreMulti+b)
plt.xlabel('Population Multipled Sim Score (dir tuning + RF location)')
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.ylabel('Population Pearson Correlation')
plt.xlim([0,1])
pearsonR, pValue = stats.pearsonr(combCombinedSimScoreMulti,combCorrMat)
plt.axis('equal')
plt.title(f'Pearson Correalation of X/Y axis = {pearsonR}, p-value = {pValue}', fontsize=8)
plt.show()


## combined loc sim score 
plt.scatter(combPairLocSimScore, combCorrMat, color='green')
# regression line
m, b = np.polyfit(combPairLocSimScore, combCorrMat, 1)
#add linear regression line to scatterplot 
plt.plot(combPairLocSimScore, m*combPairLocSimScore+b)
plt.xlabel('Population Loc Sim Score')
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.ylabel('Population Pearson Correlation')
plt.xlim([0,1])
pearsonR, pValue = stats.pearsonr(combPairLocSimScore,combCorrMat)
plt.axis('equal')
plt.title(f'Pearson Correalation of X/Y axis = {pearsonR}, p-value = {pValue}', fontsize=8)
plt.show()


## combined dir tuning sim score 
plt.scatter(combPairSimScore, combCorrMat, color='green')
# regression line
m, b = np.polyfit(combPairSimScore, combCorrMat, 1)
#add linear regression line to scatterplot 
plt.plot(combPairSimScore, m*combPairSimScore+b)
plt.xlabel('Population Dir Sim Score')
plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
plt.ylabel('Population Pearson Correlation')
plt.xlim([0,1])
pearsonR, pValue = stats.pearsonr(combPairSimScore,combCorrMat)
plt.axis('equal')
plt.title(f'Pearson Correalation of X/Y axis = {pearsonR}, p-value = {pValue}', fontsize=8)
plt.show()
