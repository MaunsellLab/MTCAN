/*
 *  OK.h
 *  OKernel
 *
 *  Copyright (c) 2018-2021 All rights reserved.
 *
 */

#import <AppKit/AppKit.h>
#import <Lablib/Lablib.h>

#define kMaxVisStim     12
#define kMaxOptoStim    2                           // NB: Index 0 is no opto stim; 1 is opto stim

typedef NS_ENUM(long, StimType) {kGaborType, kLEDType, kDelayType, kStimTypes};

enum {kLinear = 0, kLogarithmic};
//enum {kLeftArrowKey = 123, kRightArrowKey, kDownArrowKey, kUpArrowKey};

typedef struct TrialDesc {                          // description of one trial
    long    serialNumber;
    long    visualStimType;
    long    visualStimIndex;
    float   visualStimValue;
    long    visualDurMS;
    long    visualRampDurMS;
    long    gaborStartFrame;                        // start of gabor/ramp
    long    gaborRampEndFrame;                      // end of ramp
    long    gaborEndFrame;                          // end of gabor
    long    optoIndex;                              // index of delay value
    float   pulseContrast;                          // contrast of modulation [0:1]
    long    pulseDurMS;                             // duration of individual optogenetic pulses
    long    delayLEDFromTriggerMS;                  // delay after NIDAQ trigger
    float   meanPowerMW;                            // mean optogenetic power
    long    preStimMS;                              // time before stimulus onset
    BOOL    syntheticData;                          // synthetic digital input (fake mouse)
    BOOL    syntheticDigitalOutput;                 // synthetic digital output (no lab hardware)
} TrialDesc;

typedef struct BlockStatus {
    long    visualStimType;                                 // type of visual stimulus (Gabor/LED)
    long    numVisualStim;                                  // number of visual stimulus levels active
    float   stimValues[kMaxVisStim];                        // list of the contrast values
    long    numRepsEachValue[kMaxVisStim][kMaxOptoStim];    // num reps of each change
    long    repsRemainingEachValue[kMaxVisStim][kMaxOptoStim]; // num changes remaining;
    long    repsDoneThisBlock;                              // num trials done in this block
    long    repsPerBlock;                                   // num of reps in a block
    long    blocksDone;                                     // number of blocks completed
    long    blockLimit;                                     // number of blocks before stopping
} BlockStatus;

// parameters set in the Stimulus dialog

typedef struct StimSetting {
    long    stimDurMS;
    long    stimOnTimeMS;
    float   eccentricityDeg;
    float   polarAngleDeg;
    float   driftDirectionDeg;
    float   meanPowerMW;
    long    powers;
    long    powerScale;
    float   maxPowerMW;
    float   minPowerMW;
} StimSetting;

#define kOKSavedDate           @"OKSavedDate"

#ifndef    NoGlobals

#define kKNAO0CalibrationKey    @"KNAO0Calibration"
#define kKNAO1CalibrationKey    @"KNAO1Calibration"
#define kOKTaskStatusKey       @"OKTaskStatus"

// Behavior settings dialog

#define kOKBlockLimitKey            @"OKBlockLimit"
#define kOKConseqRewardLimit        @"OKConseqRewardLimit"
#define kOKRewardULKey              @"OKRewardUL"
#define kOKRespLimitMSKey           @"OKRespLimitMS"
#define kOKSubjectNumberKey         @"OKSubjectNumber"
#define kOKTooFastMSKey             @"OKTooFastMS"
#define kOKMinRunTimeS              @"OKMinRunTime"
#define kOKTotalRunTimeS            @"OKTotalRunTime"

// Stimulus settings dialog for visual stimulus

#define kOKContrastsKey            @"OKContrasts"
#define kOKContrastScaleKey        @"OKContrastScale"
#define kOKIdleGrayKey             @"OKIdleGray"
#define kOKGaborArrayKey           @"OKGaborArray"
#define kOKLEDArrayKey             @"OKLEDArray"
#define kOKMaxContrastPCKey        @"OKMaxContrastPC"
#define kOKMinContrastPCKey        @"OKMinContrastPC"
#define kOKMaxPreStimMSKey         @"OKMaxPreStimMS"
#define kOKMinPreStimMSKey         @"OKMinPreStimMS"
#define kOKPowerScaleKey           @"OKPowerScale"
#define kOKPowersKey               @"OKPowers"
#define kOKMaxPowerMWKey           @"OKMaxPowerMW"
#define kOKMinPowerMWKey           @"OKMinPowerMW"
#define kOKRampDurMSKey            @"OKRampDurMS"
#define kOKStimDurMSKey            @"OKStimDurMS"
#define kOKVisualStimTypeKey       @"OKVisualStimType"

// Stimulus settings dialog for optical stimulus

#define kOKMeanPowerMWKey       @"OKMeanPowerMW"
#define kOKPulseContrastKey     @"OKPulseContrast"
#define kOKPulseDurMSKey        @"OKPulseDurMS"
#define kOKOpticalDelayKey      @"OKOpticalDelay"
#define kOKOpticalRepsKey       @"OKOpticalReps"

// Rig hardware configuration parameters

#import "OKernel.h"

extern OKernel   *task;
extern BlockStatus  blockStatus;

#endif  /* NoGlobals */

