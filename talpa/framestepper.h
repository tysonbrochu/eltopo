
// ---------------------------------------------------------
//
//  framestepper.h
//  Tyson Brochu 2011
//  
//  Control frame sub-steps when the sim dt is not equal to the frame dt.
//
// ---------------------------------------------------------


#ifndef EL_TOPO_FRAMESTEPPER_H
#define EL_TOPO_FRAMESTEPPER_H

class FrameStepper 
{
    
    
public:
    
    double frame_length;    // length of a frame in seconds
    int frame_count;              // always the frame currently being processed
    double current_time;          // current time within a frame
    int step_count;               // current step within the frame
    
    
    FrameStepper ( double frame_len ) :
    frame_length(frame_len),
    frame_count(0),
    current_time(0), 
    step_count(0)
    {}
    
    // adjust the timestep to land on a frame time, or to use more evenly spaced steps if close to a frame time
    double get_step_length(double max_step)
    {
        if(current_time + max_step > frame_length)
            max_step = frame_length - current_time;
        
        return max_step;     
    }
    
    // we're done when current time is very close or past the frame_length
    bool done_frame()
    {
        return current_time >= frame_length - 1e-7;
    }
    
    void advance_step(double step_length)
    {
        current_time += step_length;
        ++step_count;
    }
    
    void next_frame()
    {
        current_time = 0;
        step_count = 0;
        ++frame_count;
    }   
    
    int get_step_count()
    {
        return step_count;
    }
    
    int get_frame()
    {
        return frame_count;
    }
    
    double get_time()
    {
        return (frame_count-1)*frame_length + current_time;
    }
    
};


#endif

