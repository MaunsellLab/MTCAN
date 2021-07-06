//
//  Fixate.h
//  Fixate
//
//  Created by John Maunsell on 12/23/04.
//  Copyright 2004-2021 All rights reserved.
//

#ifndef Fixate_h
#define Fixate_h

#import "FTStimuli.h"
#import "FTStateSystem.h"
#import "FTEyeXYController.h"

@interface Fixate : LLTaskPlugIn {

    NSPoint                 currentEyesUnits[kEyes];
    FTEyeXYController       *eyeXYController;                // Eye position display
    NSWindowController      *summaryController;
    NSWindowController      *xtController;
    NSArray                 *topLevelObjects;
}

@property (strong) FTStimuli *stimuli;
@property (strong) LLScheduleController *scheduler;

- (void)dataCollect:(NSTimer *)timer;
- (IBAction)doFixSettings:(id)sender;
- (IBAction)doJuice:(id)sender;
- (IBAction)doReset:(id)sender;

#ifndef    NoGlobals
extern Fixate   *task;
#endif /* NoGlobals */

@end

#endif /* Fixate_h */
