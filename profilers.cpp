#include "multi_filament_header.hpp"
#include "profilers.hpp"
#include "print_functions.hpp"

void timer::start(){
	this->start_time = std::chrono::high_resolution_clock::now();
};


int timer::end(){
	auto elapsed_time = (std::chrono::high_resolution_clock::now() - this->start_time);
	int elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time).count();
	//int elapsed_sec = elapsed_ms/1000;
	this->value += elapsed_ms;
	this->count++;
	return elapsed_ms;
};

float timer::average(){
	if (this->count == 0) {
		return 0.0;
	} else {
		return (float)this->value / (float)this->count;
	}
};

void timer::reset(){
	this->value = 0;
	this->count = 0;
};

void timer::print_total(){
	cout << "Total time " << this->name << ": "
		 << format_time(this->value/1000.) << endl;
}

string timer::print_total_short(){
	return format_time(this->value/1000.);
}

void timer::print_average(){
	cout << "Average time " << this->name << ":            "
		 << format_time(this->average()/1000.) << endl;
}

void profiler::print_all(int nt){
	cout << "[" << this->fcm_make_zero.print_total_short() << "] ";
	cout << "[" << this->fcm_gaussian_setup.print_total_short() << "] ";
	cout << "[" << this->fcm_force_distribution.print_total_short() << "] ";
	cout << "[" << this->fcm_fft_forward.print_total_short() << "] ";
	cout << "[" << this->fcm_flow_solve.print_total_short() << "] ";
	cout << "[" << this->fcm_fft_backward.print_total_short() << "] ";
	cout << "[" << this->fcm_particle_velocities_rotations.print_total_short() << "] ";
	cout << "[" << this->totalTimeFCM.print_total_short() << "] ";
	cout << "[" << this->totalTimeCopy.print_total_short() << "] ";
	cout << "[" << this->totalTimeBuildJacobian.print_total_short() << "] ";
	cout << "[" << this->totalTimeCollisionBarrier.print_total_short() << "] ";
	cout << "[" << this->totalTimeErrorCheck.print_total_short() << "] ";
	cout << "[" << this->totalTimeSolveJacobian.print_total_short() << "] ";
	cout << "[" << this->totalTimestepTime.print_total_short() << "] ";
	cout << "[[" << "\033[1;33m"
		 << "<" << format_time(this->totalTime.average()/1000.*(TimeSteps-nt))
		 << "\033[0m" << "]] ";

}

void profiler::print_names(){
	cout << "[" << this->fcm_make_zero.short_name << "] ";
	cout << "[" << this->fcm_gaussian_setup.short_name << "] ";
	cout << "[" << this->fcm_force_distribution.short_name << "] ";
	cout << "[" << this->fcm_fft_forward.short_name << "] ";
	cout << "[" << this->fcm_flow_solve.short_name << "] ";
	cout << "[" << this->fcm_fft_backward.short_name << "] ";
	cout << "[" << this->fcm_particle_velocities_rotations.short_name << "] ";
	cout << "[" << this->totalTimeFCM.short_name << "] ";
	cout << "[" << this->totalTimeCopy.short_name << "] ";
	cout << "[" << this->totalTimeBuildJacobian.short_name << "] ";
	cout << "[" << this->totalTimeCollisionBarrier.short_name << "] ";
	cout << "[" << this->totalTimeErrorCheck.short_name << "] ";
	cout << "[" << this->totalTimeSolveJacobian.short_name << "] ";
	cout << "[" << this->totalTimestepTime.short_name << "] ";
	cout << "[[" << this->totalTime.short_name << "]] ";
}
