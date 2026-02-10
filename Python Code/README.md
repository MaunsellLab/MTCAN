Arrays and Explanations

prefNormalizedWithBase: 
Neuronal responses to preferred stimuli at the 5 testing sites normalized to the preferred center response while retaining
baseline activity. The rows correspond to measurements from a single neuron and columns are ordered from center to periphery 
(first column=center, last column=offset location 4). The shape of the array is 90x5 (90 neurons, 5 conditions)
	
•	

nonprefNormalizedWithBase:
Same as above but to the non-preferred stimulus

•

pnNormalizedWithBase: same as above but responses to paired stimulus condition. Here a preferred stimulus is presented
at the center, while a non-preferred stimulus appears at 1 of 4 offset location sites. Thus in this array, the first column
corresponds to the paired condition in which there is a preferred stimulus at the center with a non-preferred stimulus
at offset location 1. Shape of this matrix is 90x4 (90 neurons, 4 paired conditions)

•	

npNormalizedWithBase: same as above but with non-preferred stimulus in the center and a preferred stimulus in the offset 
sites

•	

ppNormalizedWithBase: same as above but with preferred + preferred condition. 

•	

nnNormalizedWithBase: same as above but with non-preferred + non-preferred condition

•	

sponNormalizedToPrefCenter: spontaneous rate for each neuron normalized to the response to a preferred stimulus alone 
in the center

•	

normValBase: The normalization values (preferred center responses) that were used as divisors.
•	

transectNormalized: Transect response matrix normalized to the same preferred center. These are the full mapping stimuli
using only the preferred stimulus in the opposite section of the RF. Shape 90x5 (90 neurons, 5 conditions). The first column
is offset location -4 while last column is center (location 0). 

Notes
	•	All normalizations were performed by dividing each response matrix by the corresponding preferred center response for that session/unit.
	•	Baseline activity (spontaneous responses) was preserved in all normalizations.
	•	Stimulus naming convention:
	•	Pref = preferred direction of motion
	•	Non-pref = opposite/non-preferred direction
	•	Center / Periphery = stimulus location within the receptive field