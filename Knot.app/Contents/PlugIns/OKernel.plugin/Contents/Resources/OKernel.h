//
//  OKernel.h
//  OKernel
//
//  Copyright 2006-2021 All rights reserved.
//

#import "OKMatlabController.h"
#import "OKStateSystem.h"
#import "OKStimuli.h"

@interface OKernel:LLTaskPlugIn {

    NSWindowController          *behaviorController;
    NSWindowController          *summaryController;
    NSArray                     *topLevelObjects;

    IBOutlet NSTableColumn      *delay0Column;
    IBOutlet NSTableColumn      *delay1Column;
    IBOutlet NSTextFieldCell    *maxOpticalPowerText;
    IBOutlet NSTextFieldCell    *maxVisualPowerText;
    IBOutlet NSTabView          *visualStimTabView;
    IBOutlet NSWindow           *stimulusDialog;
}

@property (readonly, strong) LLDigitalOut *digitalOut;
@property (assign) id nidaqTask;
@property (strong) NSArray *observerKeyArray;
@property (strong) OKMatlabController *matlabController;
@property (assign) BOOL resetFlag;
@property (strong) OKStimuli *stimuli;
@property (assign) long trialCounter;

- (IBAction)doGaborSettings:(id)sender;
- (IBAction)doJuice:(id)sender;
- (IBAction)doReset:(id)sender;
- (IBAction)doSettings:(id)sender;

@end

