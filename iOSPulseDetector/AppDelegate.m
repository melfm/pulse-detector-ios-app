/*****************************************************************************
 *   AppDelegate.m
 ******************************************************************************
 *   by Melissa Mozifian, 25th Dec 2014
 ******************************************************************************
 *
 *
 *
 *
 *****************************************************************************/

#import "AppDelegate.h"

@implementation AppDelegate

@synthesize  alert;
@synthesize allowCamAccess;
@synthesize cancelAccess;

- (BOOL)application:(UIApplication *)application didFinishLaunchingWithOptions:(NSDictionary *)launchOptions
{
    // Override point for customization after application launch.
    return YES;
}
							
- (void)applicationWillResignActive:(UIApplication *)application
{
    // Sent when the application is about to move from active to inactive state. This can occur for certain types of temporary interruptions (such as an incoming phone call or SMS message) or when the user quits the application and it begins the transition to the background state.
    // Use this method to pause ongoing tasks, disable timers, and throttle down OpenGL ES frame rates. Games should use this method to pause the game.
}

- (void)applicationDidEnterBackground:(UIApplication *)application
{
    // Use this method to release shared resources, save user data, invalidate timers, and store enough application state information to restore your application to its current state in case it is terminated later. 
    // If your application supports background execution, this method is called instead of applicationWillTerminate: when the user quits.
}

- (void)applicationWillEnterForeground:(UIApplication *)application
{
    // Called as part of the transition from the background to the inactive state; here you can undo many of the changes made on entering the background.
}

- (void)applicationDidBecomeActive:(UIApplication *)application
{
    // Restart any tasks that were paused (or not yet started) while the application was inactive. If the application was previously in the background, optionally refresh the user interface.
    alert=   [UIAlertController alertControllerWithTitle:@"\"PulseDetector\" Would Like to Access the Camera"
                                                 message:nil
                                          preferredStyle:UIAlertControllerStyleAlert];
    
    allowCamAccess = [UIAlertAction
                      actionWithTitle:@"Ok"
                      style:UIAlertActionStyleDefault
                      handler:^(UIAlertAction * action)
                      {
                          [alert dismissViewControllerAnimated:YES completion:nil];

                          
                      }];
    cancelAccess = [UIAlertAction
                    actionWithTitle:@"Don't Allow"
                    style:UIAlertActionStyleDefault
                    handler:^(UIAlertAction * action)
                    {
                        [alert dismissViewControllerAnimated:YES completion:nil];
                        // Exit the app? Cancel? What to do?
                        
                    }];
    [alert addAction:cancelAccess];
    [alert addAction:allowCamAccess];
    [alert setModalPresentationStyle:UIModalPresentationPopover];

    [self.window.rootViewController presentViewController:alert animated:YES completion:nil];
    //[self presentViewController:alert animated:YES completion:nil];
}

- (void)applicationWillTerminate:(UIApplication *)application
{
    // Called when the application is about to terminate. Save data if appropriate. See also applicationDidEnterBackground:.
}

@end
