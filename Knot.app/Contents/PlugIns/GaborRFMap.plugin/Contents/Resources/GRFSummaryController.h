//
//  GRFSummaryController.h
//  Experiment
//
//  Copyright (c) 2006-2021 All rights reserved.
//

#import "GRF.h"

@interface GRFSummaryController : LLScrollZoomWindow {

    double				accumulatedRunTimeS;
	BlockStatus			blockStatus;
	long				dayComputer;			// Count of trials with computer certification errors
	long				dayEOTs[kEOTTypes];
    long				dayEOTTotal;
    NSDictionary		*fontAttr;
    NSDictionary		*labelFontAttr;
    NSDictionary		*leftFontAttr;
    double				lastStartTimeS;
	MappingBlockStatus	mappingBlockStatus;
    BOOL				newTrial;
	long				recentComputer;			// Count of trials with computer certification errors
    long				recentEOTs[kEOTTypes];
    long				recentEOTTotal;
	StimParams			stimParams;
    long 				taskMode;
//	TrialDesc			trial;

    IBOutlet			LLEOTView *dayPlot;
    IBOutlet			LLEOTHistoryView *eotHistory;
    IBOutlet			NSTableView *percentTable;
    IBOutlet			LLEOTView *recentPlot;
    IBOutlet			NSTableView *trialTable;
}

- (NSDictionary *)makeAttributesForFont:(NSFont *)font alignment:(NSTextAlignment)align tailIndex:(float)indent;
- (id)percentTableColumn:(NSTableColumn *)tableColumn row:(long)row;
- (id)trialTableColumn:(NSTableColumn *)tableColumn row:(long)row;

@end
