#ifndef _UTILITY_H
#define	_UTILITY_H

namespace data_manage{

class Utility{
	public:
		static void decrement(int& i){i--;}
		static void increment(int& i){i++;}
		
	private:
		// private so no instance can be created
		Utility();
};
}

#endif