#include "EPICLib/Output_tee_globals.h"
#include "HZ13_v5_device.h"


// for use in non-dynamically loaded models
Device_base * create_concrete_device()
{
	return new VisSearch_device("HZ13 device", Normal_out);
}

Device_base * create_HZ13_device()
{
	return new VisSearch_device("HZ13 device", Normal_out);
}

// the class factory functions to be accessed with dlsym
extern "C" Device_base * create_device() 
{
    return create_HZ13_device();
}

extern "C" void destroy_device(Device_base * p) 
{
    delete p;
}
