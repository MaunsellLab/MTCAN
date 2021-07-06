/*
 *  VF.h
 *  VFMap
 *
 *  Copyright (c) 2006-2021 All rights reserved.
 *
 */

#import <AppKit/AppKit.h>
#import <Lablib/Lablib.h>

@class VFDigitalOut;

#define kPI          		(atan(1) * 4)
#define k2PI         		(atan(1) * 4 * 2)
#define kRadiansPerDeg      (kPI / 180.0)
#define kDegPerRadian		(180.0 / kPI)

typedef struct StimDesc {
	long	sequenceIndex;
	long	stimOnFrame;
	long	stimOffFrame;
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
	long	intertrialMS;
	long	rewardSchedule;
	long	rewardMS;
} BehaviorSetting;

#ifndef	NoGlobals

// Behavior settings dialog

#define kVFDoSoundsKey             @"VFDoSounds"
#define kVFIntertrialMSKey         @"VFIntertrialMS"
#define kVFRewardMSKey             @"VFRewardMS"
#define kVFMinRewardMSKey          @"VFMinRewardMS"
#define kVFRewardScheduleKey       @"VFRewardSchedule"
#define kVFTaskStatusKey           @"VFTaskStatus"

// Stimulus Parameters

#define kVFDoDigitalKey             @"VFDoDigital"
#define kVFMapInterstimDurationMSKey @"VFMapInterstimDurationMS"
#define kVFMappingBlocksKey         @"VFMappingBlocks"
#define kVFMapStimDurationMSKey     @"VFMapStimDurationMS"
#define kVFMapTemporalModulationKey @"VFMapTemporalModulation"
#define kVFStimTablesKey            @"VFStimTables"
#define kVFStimTableCountsKey       @"VFStimTableCounts"
#define kVFMapStimContrastPCKey     @"VFMapStimContrastPC"
#define kVFMapStimRadiusSigmaRatioKey @"VFMapStimRadiusSigmaRatio"

#import "VFStimuli.h"
#import "VFMap.h"

extern BehaviorSetting				behaviorSetting;
extern MappingBlockStatus			mappingBlockStatus;
extern BOOL							resetFlag;
extern VFDigitalOut                 *digitalOut;
extern VFMap                        *task;
extern long                         trialCounter;

#endif  /* NoGlobals */
