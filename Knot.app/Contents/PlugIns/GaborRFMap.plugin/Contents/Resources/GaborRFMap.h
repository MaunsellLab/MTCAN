//
//  GaborRFMap.h
//  GaborRFMap
//
//  Copyright 2006-2021 All rights reserved.
//

#import "GRF.h"
#import "GRFStateSystem.h"
#import "GRFEyeXYController.h"
#import "GRFRoundToStimCycle.h"
#import "GRFStateSystem.h"

@class GRFMapStimTable;

@interface GaborRFMap:LLTaskPlugIn {

    NSWindowController      *behaviorController;
    NSPoint                 currentEyesUnits[kEyes];
    GRFEyeXYController      *eyeXYController;                // Eye position display
    NSWindowController      *spikeController;
    NSWindowController      *summaryController;
    NSArray                 *topLevelObjects;
}

@property (readonly, strong) LLDigitalOut *digitalOut;
@property (readonly, strong) GRFMapStimTable *mapStimTable0;
@property (readonly, strong) GRFMapStimTable *mapStimTable1;
@property (strong, readonly) LLEyeWindow *respWindow;
@property (strong) GRFStimuli *stimuli;
@property TrialDesc         trial;

- (IBAction)doFixSettings:(id)sender;
- (IBAction)doJuice:(id)sender;
- (void)doJuiceOff;
- (IBAction)doReset:(id)sender;
- (IBAction)doTaskGaborSettings:(id)sender;


- (void)updateChangeTable;

@end

