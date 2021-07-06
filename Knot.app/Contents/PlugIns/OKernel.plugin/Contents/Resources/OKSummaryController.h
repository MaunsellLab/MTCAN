//
//  OKSummaryController.h
//  Experiment
//
//  Copyright (c) 2018-2021 All rights reserved.
//

#import "OK.h"

@interface OKSummaryController : LLScrollZoomWindow {

    double				accumulatedRunTimeS;
	BlockStatus			blockStatus;
	long				dayComputer;			// Count of trials with computer certification errors
	long				dayEOTs[kEOTTypes];
    long				dayEOTTotal;
    NSDictionary		*fontAttr;
    NSDictionary		*labelFontAttr;
    NSDictionary		*leftFontAttr;
	long 				lastEOTCode;
    double				lastStartTimeS;
    long				recentEOTs[kEOTTypes];
    long				recentEOTTotal;
    BOOL				newTrial;
	long				recentComputer;			// Count of trials with computer certification errors
    long                subjectNumber;
    long 				taskMode;
	TrialDesc			trial;

    IBOutlet			LLEOTView *dayPlot;
    IBOutlet			LLEOTHistoryView *eotHistory;
    IBOutlet			NSTableView *percentTable;
    IBOutlet			LLEOTView *recentPlot;
    IBOutlet			LLSymbolView *symbolView;
    IBOutlet			NSTableView *trialTable;
}

- (NSDictionary *)makeAttributesForFont:(NSFont *)font alignment:(NSTextAlignment)align tailIndex:(float)indent;
- (id)percentTableColumn:(NSTableColumn *)tableColumn row:(long)row;
- (id)trialTableColumn:(NSTableColumn *)tableColumn row:(long)row;

@end
