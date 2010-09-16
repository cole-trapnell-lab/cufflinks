#ifndef PROGRESS_H
#define PROGRESS_H

#include "time.h"
using namespace std;

const int BAR_BUF_SIZE = 28;
const char SYMBOL = '*';

class ProgressBar
{
	char _bar_buf[BAR_BUF_SIZE];
	string _process;
	int _num_complete;
	int _tot_num;
public:
	ProgressBar(string process, int tot_num) 
	{ 
		_tot_num = tot_num;
		_process = process;
		_num_complete = -1;
		for(int i=0; i < BAR_BUF_SIZE; ++i) _bar_buf[i] = ' ';
		_bar_buf[0] = '[';
		_bar_buf[BAR_BUF_SIZE-2] = ']';
		_bar_buf[BAR_BUF_SIZE-1] = '\0';
		
		
		time_t rawtime;
  		struct tm * timeinfo;
  		char time_buf [80];

  		time ( &rawtime );
  		timeinfo = localtime ( &rawtime );
		
		strftime (time_buf,80,"%H:%M:%S",timeinfo);
		fprintf(stderr, "[%s] %s\n", time_buf, _process.c_str());
	}
	
	void update(char* bundle_label_buf, int inc_amt)
	{
		_num_complete += inc_amt;
		int percent = (_num_complete * 100)/_tot_num;	
		int last_bar = percent/(100/(BAR_BUF_SIZE-3));
		for (int i=1; i <= last_bar; ++i)
			_bar_buf[i] = SYMBOL;
		
		char line_buf[81];
		bundle_label_buf[30] = '\0';
		sprintf(line_buf, "\rProcessing Locus %-30s %s %3d%%", bundle_label_buf, _bar_buf, percent);
		fprintf(stderr,"%s",line_buf);
	}
	
	void complete()
	{
		for (int i=1; i < BAR_BUF_SIZE-2; ++i)
			_bar_buf[i] = SYMBOL;
		fprintf(stderr, "\r%-47s %s %3d%%\n", "Completed", _bar_buf, 100);
	}
};

#endif
