//
//  GRFBehaviorController.h
//  GaborRFMap
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "GRF.h"

#define kEOTs				(kEOTIgnored + 2)				// Plot EOT types up to ignored, + 1 for probes
#define kEOTProbeCorrect	(kEOTIgnored + 1)
#define kMaxRT				1000

@interface GRFBehaviorController : LLScrollZoomWindow {

	NSView			*documentView;
    LLHistView		*hist[kMaxOriChanges];
    LLViewScale		*histScaling;
    NSColor 		*highlightColor;
    long			histHighlightIndex;
    NSMutableArray	*labelArray;
	BlockStatus		lastBlockStatus;
	long			lastMaxTargetTimeMS;
	long			lastMinTargetTimeMS;
	long			maxTargetTimeMS;
	long			minTargetTimeMS;
    NSMutableArray	*performance[kEOTs];				// an array of LLBinomDist
    NSMutableArray	*performanceByTime[kEOTs];			// an array of LLBinomDist
    LLPlotView		*perfPlot;
    LLPlotView		*perfTimePlot;
    LLPlotView		*reactPlot;
    NSMutableArray	*reactTimes;						// an array of LLNormDist
	long			saccadeStartTimeMS;
    long			responseTimeMS;
    double			rtDist[kMaxOriChanges][kMaxRT];
	long            stimStartTimeMS;
	long            targetOnTimeMS;
    NSMutableArray	*timeLabelArray;
    NSMutableArray	*xAxisLabelArray;
}

- (void)changeResponseTimeMS;
- (void)checkTimeParams;
- (void)checkParams;
- (LLHistView *)myInitHist:(LLViewScale *)scale data:(double *)data;
- (void)makeLabels;
- (void)makeTimeLabels;
- (void)positionPlots;
- (void)reset:(NSData *)eventData eventTime:(NSNumber *)eventTime;

@end
