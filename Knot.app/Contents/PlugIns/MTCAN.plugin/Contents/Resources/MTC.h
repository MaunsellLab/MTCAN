/*
 *  MTC.h
 *  MTCAN
 *
 *  Copyright (c) 2006-2021 All rights reserved.
 *
 */

#import <AppKit/AppKit.h>
#import <Lablib/Lablib.h>


// The following should be changed to be unique for each application

//typedef enum {kFirstGabor = 0, kTaskGabor = 0, kMapGabor0, kMapGabor1, kGabors} GaborType;
//typedef enum {kFirstGabor = 0, kTaskGabor = 0, kGabors} GaborType;

typedef NS_ENUM (long, LocType) {kGabor0 = 0, kGabor1, kGaborFar0, kGaborFar1, kNumLocations};
typedef NS_ENUM (short, ListType)  {kNilStim = 0, kValidStim, kTargetStim, kFrontPadding, kBackPadding, kListTypes};
typedef NS_ENUM (short, RFStimType) {kNoStim = 0, kNullStim, kPrefStim, kRFStimTypes};

enum {kLinear = 0, kLogarithmic};
enum {kUniform = 0, kExponential};
enum {kAuto = 0, kManual};
enum {kRewardFixed = 0, kRewardVariable};
enum {kMyEOTCorrect = 0, kMyEOTMissed, kMyEOTEarlyToValid, kMyEOTEarlyToInvalid, kMyEOTBroke,
				kMyEOTIgnored, kMyEOTQuit, kMyEOTTypes};
enum {
//        kTrialStartDigitOutCode =  0x0010,
//        kFixateDigitOutCode =      0x0020,
//        kStimulusOnDigitOutCode =  0x0030,
//        kStimulusOffDigitOutCode = 0x0040,
        kTargetOnDigitOutCode =    0x0050,
//        kSaccadeDigitOutCode =     0x0060,
//        kTrialEndDigitOutCode =    0x0070,
//        kMTCANDigitOutCode =  0x0080,
        kFixOnDigitOutCode =       0x0090};

#define	kMaxDirChanges      12
#define kPI                 (atan(1) * 4)
#define kRadiansPerDeg      (kPI / 180.0)
#define kDegPerRadian       (180.0 / kPI)
#define kStimPerLocation    (kRFStimTypes * kRFStimTypes)

typedef struct {
	long	levels;				// number of active stimulus levels
	float   maxValue;			// maximum stimulus value (i.e., direction change in degree)
	float   minValue;			// minimum stimulus value
} StimParams;

typedef struct StimDesc {
    ListType listTypes[kNumLocations];        // list type (valid, target, padding, etc.)
    RFStimType stimTypes[kNumLocations];      // stim type (pref, null, or intermediate)
	long	stimOnFrame;
	long	stimOffFrame;
    long    frameOnRendered;
    long    frameOffRendered;
	float	azimuthDeg[kNumLocations];
	float	elevationDeg[kNumLocations];
    float   directionsDeg[kNumLocations];   // direction at each location
    float   contrasts[kNumLocations];       // contrast value at each location
} StimDesc;

typedef struct TrialDesc {
	BOOL	instructTrial;
	BOOL	catchTrial;
    LocType attendLoc;                      // location cued for attention
    LocType targetLoc;                      // location of the target
    long    targetIndex;                    // position of target in stimulus sequence
	long	targetOnTimeMS;				    // time from first stimulus (start of stimlist) to the target
    long    numStim;                        // number of stimuli in stimList
	long	directionChangeIndex;
	float	directionChangeDeg;
} TrialDesc;

typedef struct BlockStatus {
	long	changes;
	float	directionChangeDeg[kMaxDirChanges];
	long	validReps[kMaxDirChanges];
	long	validRepsDone[kMaxDirChanges];
	long	invalidReps[kMaxDirChanges];
	long	invalidRepsDone[kMaxDirChanges];
	long	instructTrials;			        // number of instruction trials to be done
    long    stimPerLocation;                // reps of each stimulus at one location before changing location
    long    blockLimit;
    long    attendLoc;                      // currently cued site
    long    instructDone;                   // number of instruction trials done
    long    stimDone[kNumLocations];        // number of stimuli done in the current block
    long    locationsDone;                  // number of blocks to complete at one location before moving to another
    long    blocksDone;
} BlockStatus;

// put parameters set in the behavior controller
typedef struct BehaviorSetting {
	long	blocks;
	long	intertrialMS;
	long	acquireMS;
	long	fixGraceMS;
	long	fixateMS;
	long	fixateJitterPC;
	long	responseTimeMS;
	long	tooFastMS;
	long	minSaccadeDurMS;
	long	breakPunishMS;
	long	rewardMS;
	float	fixWinWidthDeg;
	float	respWinWidthDeg;
} BehaviorSetting;

// put parameters set in the Stimulus controller
typedef struct StimSetting {
	long	stimDurationMS;
	long	stimDurJitterPC;
	long	interStimMS;
	long	interStimJitterPC;
	long	stimLeadMS;
	float	stimSpeedHz;
	long	stimDistribution;
	long	minTargetOnTimeMS;
	long	meanTargetOnTimeMS;
	long	maxTargetOnTimeMS;
	float	eccentricityDeg;
	float	polarAngleDeg;
	float	driftDirectionDeg;
	float	contrastPC;
	long	changeRemains;
} StimSetting;

#ifndef	NoGlobals

// Behavior settings dialog
#define kMTCAcquireMSKey            @"MTCAcquireMS"
#define kMTCAlphaTargetDetectionTaskKey @"MTCAlphaTargetDetectionTask"
#define kMTCBlockLimitKey           @"MTCBlockLimit"
#define kMTCCatchTrialPCKey         @"MTCCatchTrialPC"
#define kMTCCatchTrialMaxPCKey      @"MTCCatchTrialMaxPC"
#define kMTCCueMSKey                @"MTCCueMS"
#define kMTCFixateMSKey             @"MTCFixateMS"
#define kMTCFixateOnlyKey           @"MTCFixateOnly"
#define kMTCFixJitterPCKey          @"MTCFixJitterPC"
#define kMTCFixWindowWidthDegKey    @"MTCFixWindowWidthDeg"
#define kMTCInstructionTrialsKey    @"MTCInstructionTrials"
#define kMTCIntertrialMSKey         @"MTCIntertrialMS"
#define kMTCMinTargetMSKey          @"MTCMinTargetMS"
#define kMTCMaxTargetMSKey          @"MTCMaxTargetMS"
#define kMTCMeanTargetMSKey         @"MTCMeanTargetMS"
#define kMTCNontargetContrastPCKey  @"MTCNontargetContrastPC"
#define kMTCRespSpotSizeDegKey      @"MTCRespSpotSizeDeg"
#define kMTCRespWindowWidthDegKey   @"MTCRespWindowWidthDeg"
#define kMTCRandTaskGaborDirectionKey @"MTCRandTaskGaborDirection"
#define kMTCSaccadeTimeMSKey        @"MTCSaccadeTimeMS"
#define kMTCStimPerLocationKey      @"MTCStimPerLocation"
#define kMTCStimDistributionKey     @"MTCStimDistribution"
#define kMTCTaskStatusKey           @"MTCTaskStatus"

// Stimulus Parameters
#define kMTCChangePersistsKey       @"MTCChangePersists"
#define kMTCDirectionDegKey         @"MTCDirectionDeg"
#define kMTCInterstimJitterPCKey    @"MTCInterstimJitterPC"
#define kMTCInterstimJitterTauMSKey @"kMTCInterstimJitterTauMS"
#define kMTCInterstimMSKey          @"MTCInterstimMS"
#define kMTCInvalidTrialPCKey       @"MTCInvalidTrialPC"
#define kMTCMapInterstimDurationMSKey @"MTCMapInterstimDurationMS"
#define kMTCMappingBlocksKey        @"MTCMappingBlocks"
#define kMTCMapStimDurationMSKey    @"MTCMapStimDurationMS"
#define kMTCStimDurationMSKey       @"MTCStimDurationMS"
#define kMTCStimJitterPCKey         @"MTCStimJitterPC"
#define kMTCChangeRemainKey         @"MTCChangeRemain"
#define kMTCChangeArrayKey          @"MTCChangeArray"
#define kMTCRelDistContrastKey      @"MTCRelDistContrast"
#define kMTCStimLeadMSKey           @"MTCStimLeadMS"
#define kMTCStimRepsPerLocationKey  @"MTCStimRepsPerLocation"
#define kMTCStimTablesKey           @"MTCStimTables"
#define kMTCStimTableCountsKey      @"MTCStimTableCounts"
#define kMTCMapStimContrastPCKey    @"MTCMapStimContrastPC"
#define kMTCTargetAlphaKey          @"MTCTargetAlpha"
#define kMTCTargetRadiusKey         @"MTCTargetRadius"

extern long		argRand;

#import "MTCStimuli.h"
#import "MTCAN.h"

extern BehaviorSetting	    behaviorSetting;
extern BOOL				    resetFlag;
extern MTCAN                *task;
extern long                 trialCounter;

#endif  /* NoGlobals */
