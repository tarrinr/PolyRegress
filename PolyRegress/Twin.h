#ifndef __T_WIN__
#define __T_WIN__

#include <windows.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Twin {
public:
	Twin(std::string);
	int rows();
	int columns();
	template<typename T>
	void println(T);
	void println();
	template<typename T>
	void print(T);
	template<typename T>
	void printmulti(T);
	void display();
	void clear();
private:
	void size();
	void line();
	void wincheck();
	int x,y;
	std::vector<std::string> buffer;
	std::string pname;

};

void Twin::size() {

	CONSOLE_SCREEN_BUFFER_INFO csbi;

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
	x = csbi.srWindow.Right - csbi.srWindow.Left + 1;
	y = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;

}

int Twin::rows() {
	size();
	return y;
}

int Twin::columns() {
	size();
	return x;
}

template<typename T>
void Twin::print(T in) {

	if (buffer.size() == 0) buffer.push_back("");

	std::ostringstream oss;
	oss << in;

	std::string temp = buffer.back() + oss.str();

	buffer.pop_back();
	buffer.push_back(temp);

}

template<typename T>
void Twin::println(T in) {

	std::ostringstream oss;
	oss << in;

	buffer.push_back(oss.str());

}

template<typename T>
void Twin::printmulti(T in) {

	std::ostringstream oss;

	int i = 0;

	while (oss << in[i]) {
		buffer.push_back(oss.str());
		i++
	}
}

void Twin::println() {
	buffer.push_back("");
}

void Twin::display() {

	wincheck();

	line();
	std::cout << "  " << pname << std::endl;
	line();
	std::cout << std::endl;

	for (auto&& i : buffer) {
		if (i != "") std::cout << "  > " << i << std::endl;
		else std::cout << std::endl;
	}

	for (int i = 0; i < y - 7 - int(buffer.size()); i++) std::cout << std::endl;
	line();
	std::cout << "  " << "Input:" << std::endl << "  > ";
	buffer.clear();

}

void Twin::wincheck() {

	size();

	int b = 0;

	for (auto&& i : buffer) {

		if (int(i.length()) > b) b = i.length();

	}

	while (y - 8 - int(buffer.size()) < 1 || int(pname.length()) > x - 4 || int(b) > x - 8) {

		std::cout << "  Window too small. . .";

		std::cin.ignore(1000, '\n');
		std::cin.clear();
		size();

	}

}

Twin::Twin(std::string name) {

	pname = name;

}

void Twin::line() {

	size();

	std::cout << " ";
	for (int i = 0; i < x - 2; i++) std::cout << "=";
	std::cout << " ";

}

void Twin::clear() {
	buffer.clear();
}
#endif