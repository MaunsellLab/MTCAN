//
//  GRFSpikeController.h
//  GaborRFMap
//
//  Copyright (c) 2006-2012. All rights reserved.
//

#import "GRF.h"

#define kMaxSpikeMS				10000
#define	kLocations				3
#define kNumSpikeChannels       2

enum PlotTypes {kDirectionPlot = 0, kContrastPlot, kTFPlot, kSFPlot, kSigmaPlot, kNumPlots};

@interface GRFSpikeController : LLScrollZoomWindow {

//    NSMutableArray	*attRates;								// an array of LLNormDist
//	BlockStatus		blockStatus;
//	NSView			*documentView;
//    LLHeatMapView	*heatMaps[kNumSpikeChannels];
//    NSMutableArray	*heatMapRates[kNumSpikeChannels];             //  LLNormDist for plotting
//    NSColor 		*highlightColor;
//    long			interstimDurMS;
//    NSMutableArray	*labelArray;
//	BlockStatus		lastBlockStatus;
//	StimParams		lastStimParams;
//    MapSettings     mapSettings[kNumSpikeChannels];
//    NSMutableArray	*rates[kNumSpikeChannels][kNumPlots];          //  LLNormDist for plotting
//    LLPlotView		*ratePlots[kNumSpikeChannels][kNumPlots];
//	unsigned		referenceOnTimeMS;                      // onset time of reference direction
//    NSMutableArray  *stimDescs[kNumSpikeChannels];
////    NSMutableArray  *stimTimes[kNumSpikeChannels];
//    long			stimDurMS;
//	float			spikePeriodMS;
//	unsigned long	targetOnTimeMS;
//	long            trialStartTime;
//	short           trialSpikes[kNumSpikeChannels][kMaxSpikeMS];
//    NSMutableArray	*xAxisLabels[kNumSpikeChannels][kNumPlots];
//    NSMutableArray	*xHMAxisLabels[kNumSpikeChannels];
//    NSMutableArray	*yHMAxisLabels[kNumSpikeChannels];
}

- (void)checkParams;
- (void)reset:(NSData *)eventData eventTime:(NSNumber *)eventTime;

@end

