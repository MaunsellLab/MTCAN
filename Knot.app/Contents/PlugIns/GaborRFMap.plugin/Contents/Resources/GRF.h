/*
 *  GRF.h
 *  GaborRFMap
 *
 *  Copyright (c) 2006-2021 All rights reserved.
 *
 */

#import <AppKit/AppKit.h>
#import <Lablib/Lablib.h>

#define kPI          		(atan(1) * 4)
#define k2PI         		(atan(1) * 4 * 2)
#define kRadiansPerDeg      (kPI / 180.0)
#define kDegPerRadian		(180.0 / kPI)

// The following should be changed to be unique for each application

typedef enum {kFirstGabor = 0, kTaskGabor = 0, kMapGabor0, kMapGabor1, kGabors} GaborType;

enum {kAttend0 = 0, kAttend1, kLocations};
enum {kLinear = 0, kLogarithmic};
enum {kUniform = 0, kExponential};
enum {kAuto = 0, kManual};
enum {kRewardFixed = 0, kRewardVariable};
enum {kNullStim = 0, kValidStim, kTargetStim, kFrontPadding, kBackPadding};
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
//        kGaborRFMapDigitOutCode =  0x0080,
        kFixOnDigitOutCode =       0x0090};

#define	kMaxOriChanges      12

typedef struct {
	long	levels;				// number of active stimulus levels
	float   maxValue;			// maximum stimulus value (i.e., direction change in degree)
	float   minValue;			// minimum stimulus value
} StimParams;

typedef struct StimDesc {
    short   stimType;
	long	gaborIndex;
	long	sequenceIndex;
	long	stimOnFrame;
	long	stimOffFrame;
    long    frameRendered;
    long    frameRateHz;
	float	orientationChangeDeg;
	float	contrastPC;
	float	azimuthDeg;
	float	elevationDeg;
	float	sigmaDeg;
	float	radiusDeg;
	float	spatialFreqCPD;
	float	directionDeg;
    float   temporalFreqHz;
	long	azimuthIndex;
	long	elevationIndex;
	long	sigmaIndex;
	long	spatialFreqIndex;
	long	directionIndex;
	long	contrastIndex;
    long    temporalFreqIndex;
    long    temporalModulation;
} StimDesc;

typedef struct TrialDesc {
	BOOL	instructTrial;
	BOOL	catchTrial;
	long	numTaskStim;
	long	targetIndex;				// index (count) of target in stimulus sequence
	long	targetOnTimeMS;				// time from first stimulus (start of stimlist) to the target
	long	orientationChangeIndex;
	float	orientationChangeDeg;
} TrialDesc;

typedef struct BlockStatus {
	long	changes;
	float	orientationChangeDeg[kMaxOriChanges];
	long	validReps[kMaxOriChanges];
	long	validRepsDone[kMaxOriChanges];
	long	invalidReps[kMaxOriChanges];
	long	invalidRepsDone[kMaxOriChanges];
	long	instructDone;			// number of instruction trials done
	long	instructTrials;			// number of instruction trials to be done
	long	sidesDone;				// number of sides (out of kLocations) done
	long	blockLimit;				// number of blocks before stopping
	long	blocksDone;				// number of blocks completed
} BlockStatus;

typedef struct MappingBlockStatus {
	long	stimDone;				// number of stim done in this block
	long	stimLimit;				// number of stim in block
	long	blocksDone;				// number of blocks completed
	long	blockLimit;				// number of blocks before stopping
} MappingBlockStatus;

typedef struct  MapParams {
    long    n;                      // number of different conditions
    float   minValue;               // smallest value tested
    float   maxValue;               // largest value tested
} MapParams;

typedef struct  MapSettings {
    MapParams    azimuthDeg;
    MapParams    elevationDeg;
    MapParams    directionDeg;
    MapParams    spatialFreqCPD;
    MapParams    sigmaDeg;
    MapParams    contrastPC;
    MapParams    temporalFreqHz;
} MapSettings;

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
	short	numberOfSurrounds;
	long	changeScale;
	long	orientationChanges;
	float	maxChangeDeg;
	float	minChangeDeg;
	long	changeRemains;
} StimSetting;


#ifndef	NoGlobals

// Behavior settings dialog

#define kGRFAcquireMSKey            @"GRFAcquireMS"
#define kGRFAlphaTargetDetectionTaskKey @"GRFAlphaTargetDetectionTask"
#define kGRFBlockLimitKey           @"GRFBlockLimit"
#define kGRFChangeScaleKey          @"GRFChangeScale"
#define kGRFCatchTrialPCKey         @"GRFCatchTrialPC"
#define kGRFCatchTrialMaxPCKey      @"GRFCatchTrialMaxPC"
#define kGRFCueMSKey                @"GRFCueMS"
#define kGRFFixateMSKey             @"GRFFixateMS"
#define kGRFFixateOnlyKey           @"GRFFixateOnly"
#define kGRFFixJitterPCKey          @"GRFFixJitterPC"
#define kGRFFixWindowWidthDegKey    @"GRFFixWindowWidthDeg"
#define kGRFInstructionTrialsKey    @"GRFInstructionTrials"
#define kGRFIntertrialMSKey         @"GRFIntertrialMS"
#define kGRFMinTargetMSKey          @"GRFMinTargetMS"
#define kGRFMaxTargetMSKey          @"GRFMaxTargetMS"
#define kGRFMeanTargetMSKey         @"GRFMeanTargetMS"
#define kGRFNontargetContrastPCKey  @"GRFNontargetContrastPC"
#define kGRFRespSpotSizeDegKey      @"GRFRespSpotSizeDeg"
#define kGRFRespWindowWidthDegKey   @"GRFRespWindowWidthDeg"
#define kGRFRandTaskGaborDirectionKey @"GRFRandTaskGaborDirection"
#define kGRFSaccadeTimeMSKey        @"GRFSaccadeTimeMS"
#define kGRFStimRepsPerBlockKey     @"GRFStimRepsPerBlock"
#define kGRFStimDistributionKey     @"GRFStimDistribution"
#define kGRFTaskStatusKey           @"GRFTaskStatus"

// Stimulus Parameters

#define kGRFInterstimJitterPCKey    @"GRFInterstimJitterPC"
#define kGRFInterstimMSKey          @"GRFInterstimMS"
#define kGRFMapInterstimDurationMSKey @"GRFMapInterstimDurationMS"
#define kGRFMappingBlocksKey        @"GRFMappingBlocks"
#define kGRFMapStimDurationMSKey    @"GRFMapStimDurationMS"
#define kGRFStimDurationMSKey       @"GRFStimDurationMS"
#define kGRFStimJitterPCKey         @"GRFStimJitterPC"
#define kGRFOrientationChangesKey   @"GRFOrientationChanges"
#define kGRFMaxDirChangeDegKey      @"GRFMaxDirChangeDeg"
#define kGRFMinDirChangeDegKey      @"GRFMinDirChangeDeg"
#define kGRFChangeRemainKey         @"GRFChangeRemain"
#define kGRFChangeArrayKey          @"GRFChangeArray"
#define kGRFStimTablesKey           @"GRFStimTables"
#define kGRFStimTableCountsKey      @"GRFStimTableCounts"
#define kGRFMapStimContrastPCKey    @"GRFMapStimContrastPC"
#define kGRFMapStimRadiusSigmaRatioKey @"GRFMapStimRadiusSigmaRatio"
#define kGRFTargetAlphaKey          @"GRFTargetAlpha"
#define kGRFTargetRadiusKey         @"GRFTargetRadius"

#define kGRFHideLeftKey             @"GRFHideLeft"
#define kGRFHideRightKey            @"GRFHideRight"
#define kGRFHideLeftDigitalKey      @"GRFHideLeftDigital"
#define kGRFHideRightDigitalKey     @"GRFHideRightDigital"
#define kGRFConvertToGratingKey     @"GRFConvertToGrating"
// Is there an ITC-18 dedicated to outputing digital words (to be Cerebus), or is it also used for controlling
// juice rewards, etc.
#define kGRFDedicatedITC18Key       @"GRFDedicatedITC18"

#define kGRFHideTaskGaborKey        @"GRFHideTaskGabor"
#define kGRFIncludeCatchTrialsinDoneListKey @"GRFIncludeCatchTrialsinDoneList"
#define kGRFMapTemporalModulationKey @"GRFMapTemporalModulation"
#define kGRFConvertToPlaidKey       @"GRFConvertToPlaid"

// Visual Stimulus Parameters

#define kGRFSpatialPhaseDegKey      @"GRFSpatialPhaseDeg"
#define kGRFTemporalFreqHzKey       @"GRFTemporalFreqHz"

// Keys for change array

#define kGRFChangeKey               @"change"
#define kGRFValidRepsKey            @"validReps"
#define kGRFInvalidRepsKey          @"invalidReps"

extern long		argRand;

#import "GRFStimuli.h"
#import "GaborRFMap.h"

extern BlockStatus					blockStatus;
extern BehaviorSetting				behaviorSetting;
extern MappingBlockStatus			mappingBlockStatus;
extern BOOL							resetFlag;
extern GaborRFMap                   *task;
extern long                         trialCounter;

#endif  /* NoGlobals */
