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
	long double _num_complete;
	long double _tot_num;
	int _num_updates;
	int _num_remaining;
	
public:
	ProgressBar() {}
	
	ProgressBar(string process, double tot_num) 
	{ 
		_tot_num = tot_num;
		_process = process;
		_num_complete = -1.0;
		_num_remaining = -1.0;
		_num_updates = 0;
		
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
	
	void update(const char* bundle_label_buf, double inc_amt)
	{
		
		_num_complete += inc_amt;
		_num_updates ++;
		
		if (cuff_verbose||cuff_quiet||_tot_num==0) return;

		char bundle_buf[28];
		bundle_buf[27] = '\0';
		strncpy(bundle_buf, bundle_label_buf, 27);
		
		int percent = (_num_complete * 100)/_tot_num;

		percent = min(percent, 99);
		
		int last_bar = percent/(100/(BAR_BUF_SIZE-3));
		for (int i=1; i <= last_bar; ++i)
			_bar_buf[i] = SYMBOL;
		
		char line_buf[82];
		snprintf(line_buf, 82, "\r> Processing Locus %-27s %s %3d%%", bundle_buf, _bar_buf, percent);
		
		fprintf(stderr,"%s",line_buf);
	}
	
	void remaining(int num_remaining)
	{
		if (cuff_verbose||cuff_quiet||_tot_num==0||num_remaining==_num_remaining) return;
		
		_num_remaining = num_remaining;
		
		int percent = (_num_complete * 100)/_tot_num;
		percent = min(percent, 99);
		
		char msg_buff[45];
		sprintf(msg_buff, "Waiting for %d threads to complete.", num_remaining);
		
		int last_bar = percent/(100/(BAR_BUF_SIZE-3));
		for (int i=1; i <= last_bar; ++i)
			_bar_buf[i] = SYMBOL;
		
		char line_buf[82];
		snprintf(line_buf, 81, "\r> %-44s %s %3d%%", msg_buff, _bar_buf, percent);
		
		fprintf(stderr,"%s",line_buf);
	}
	
	void complete()
	{
		for (int i=1; i < BAR_BUF_SIZE-2; ++i)
			_bar_buf[i] = SYMBOL;
		char complete_buf[45];
		snprintf(complete_buf, 44, "Processed %d loci.", _num_updates); 
		if (cuff_verbose||cuff_quiet)
			fprintf(stderr, "%-44s\n", complete_buf);
		else
			fprintf(stderr, "\r> %-44s %s %3d%%\n", complete_buf, _bar_buf, 100);
	}
};

#endif
