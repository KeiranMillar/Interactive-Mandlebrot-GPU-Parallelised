// Data Structures and Algorithms II : Intro to AMP and benchmarking exercise
// Ruth Falconer  <r.falconer@abertay.ac.uk>
// Adapted from C++ AMP book http://ampbook.codeplex.com/license.

#include <conio.h>
#include <windows.h>

#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <amp.h>
#include <fstream>
#include <time.h>
#include <string>
#include <amp_math.h>

#define SIZE 1<<20 // same as 2^20
#define TS 32

const int WIDTH = 640;
const int HEIGHT = 480;

// Need to access the concurrency libraries 
using namespace concurrency;

// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::cout;
using std::cin;
using std::endl;
using std::ofstream;

// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_serial_clock;
typedef std::chrono::steady_clock the_amp_clock;

const int MAX_ITERATIONS = 500;
uint32_t image[HEIGHT][WIDTH];

// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}


void report_accelerator(const accelerator a)
{
	const std::wstring bs[2] = { L"false", L"true" };
	std::wcout << ": " << a.description << " "
		<< endl << "       device_path                       = " << a.device_path
		<< endl << "       dedicated_memory                  = " << std::setprecision(4) << float(a.dedicated_memory) / (1024.0f * 1024.0f) << " Mb"
		<< endl << "       has_display                       = " << bs[a.has_display]
		<< endl << "       is_debug                          = " << bs[a.is_debug]
		<< endl << "       is_emulated                       = " << bs[a.is_emulated]
		<< endl << "       supports_double_precision         = " << bs[a.supports_double_precision]
		<< endl << "       supports_limited_double_precision = " << bs[a.supports_limited_double_precision]
		<< endl;
}

// List and select the accelerator to use
void list_accelerators()
{
	//get all accelerators available to us and store in a vector so we can extract details
	std::vector<accelerator> accls = accelerator::get_all();

	// iterates over all accelerators and print characteristics
	for (int i = 0; i<accls.size(); i++)
	{
		accelerator a = accls[i];
		report_accelerator(a);

	}

	accelerator acc = accelerator(accelerator::default_accelerator);
	std::wcout << " default acc = " << acc.description << endl;

} // list_accelerators

  // query if AMP accelerator exists on hardware
void query_AMP_support()
{
	std::vector<accelerator> accls = accelerator::get_all();
	if (accls.empty())
	{
		cout << "No accelerators found that are compatible with C++ AMP" << std::endl;
	}
	else
	{
		cout << "Accelerators found that are compatible with C++ AMP" << std::endl;
		list_accelerators();
	}
} // query_AMP_support

  // using our own structure as Complex function not available in the Concurrency namespace
struct Complex {
	float x;
	float y;
};

Complex c_add(Complex c1, Complex c2) restrict(cpu, amp) // restrict keyword - able to execute this function on the GPU and CPU
{
	Complex tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a + c;
	tmp.y = b + d;
	return tmp;
} // c_add

float c_abs(Complex c) restrict(cpu, amp)
{
	return concurrency::fast_math::sqrt(c.x*c.x + c.y*c.y);
} // c_abs

Complex c_mul(Complex c1, Complex c2) restrict(cpu, amp)
{
	Complex tmp;
	float a = c1.x;
	float b = c1.y;
	float c = c2.x;
	float d = c2.y;
	tmp.x = a*c - b*d;
	tmp.y = b*c + a*d;
	return tmp;
} // c_mul

void compute_mandlebrot(float left, float right, float top, float bottom)
{
	uint32_t *pimage = &(image[0][0]);
	extent<2> e(HEIGHT, WIDTH);
	concurrency::array_view<uint32_t, 2> av3(e, pimage);
	try
	{
		concurrency::parallel_for_each(av3.extent.tile<TS, TS>(), [=](concurrency::tiled_index<TS, TS> t_idx)  restrict(amp)
		{
			int x, y;
			index<2> idx = t_idx.global;
			//uint32_t colour = 0x010000;
			y = idx[0];
			x = idx[1];
			Complex c, z;

			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			c.x = (left + (x * (right - left) / WIDTH));
			c.y = (top + (y * (bottom - top) / HEIGHT));

			// Start off z at (0, 0).
			z.x = 0;
			z.y = 0;

			// Iterate z = z^2 + c until z moves more than 2 units

			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (c_abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = c_add(c_mul(z, z), c);

				++iterations;
			}

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				av3[idx] = 0x000000;
			}
			else
			{
				// z escaped within less than MAX_ITERATIONS
				// iterations. This point isn't in the set.

				av3[idx] = 0xFFFFFF;
			}
		});
		av3.synchronize();
	}
	catch (const Concurrency::runtime_exception& ex)
	{
		MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
	}
}

// Attempt at tiling
void time_to_tile(const int size, const std::vector<int>& v1, const std::vector<int>& v2, std::vector<int>& v3)
{
	concurrency::array_view<const int> av1(size, v1);
	concurrency::array_view<const int> av2(size, v2);
	extent<1> e(size);
	concurrency::array_view<int> av3(e, v3);
	av3.discard_data();

	// start clock for GPU version after array allocation
	the_amp_clock::time_point start = the_amp_clock::now();
	// It is wise to use exception handling here - AMP can fail for many reasons
	// and it useful to know why (e.g. using double precision when there is limited or no support)
	try
	{
		concurrency::parallel_for_each(av3.extent.tile <TS>(), [=](concurrency::tiled_index<TS> t_idx) restrict(amp)
		{
			index<1> idx = t_idx.global;
			av3[t_idx] = av1[t_idx] + av2[t_idx];
		});
		av3.synchronize();
	}
	catch (const Concurrency::runtime_exception& ex)
	{
		MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
	}
	// Stop timing
	the_amp_clock::time_point end = the_amp_clock::now();
	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Adding the vectors using tiling (data transfer and compute) takes " << time_taken << " ms." << endl;
} // vector_add_amp

  // Unaccelerated CPU element-wise addition of arbitrary length vectors using C++ 11 container vector class
void vector_add(const int size, const std::vector<int>& v1, const std::vector<int>& v2, std::vector<int>& v3)
{
	//start clock for serial version
	the_serial_clock::time_point start = the_serial_clock::now();
	//use iterator to loop over and add vector elements
	for (auto i : v3)
	{
		v3[i] = v2[i] + v1[i];
	}
	the_serial_clock::time_point end = the_serial_clock::now();
	//Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Adding the vectors serially " << time_taken << " ms." << endl;
} // vector_add

  // Accelerated  element-wise addition of arbitrary length vectors using C++ 11 container vector class
void vector_add_amp(const int size, const std::vector<int>& v1, const std::vector<int>& v2, std::vector<int>& v3)
{
	concurrency::array_view<const int> av1(size, v1);
	concurrency::array_view<const int> av2(size, v2);
	extent<1> e(size);
	concurrency::array_view<int> av3(e, v3);
	av3.discard_data();

	// start clock for GPU version after array allocation
	the_amp_clock::time_point start = the_amp_clock::now();
	// It is wise to use exception handling here - AMP can fail for many reasons
	// and it useful to know why (e.g. using double precision when there is limited or no support)
	try
	{
		concurrency::parallel_for_each(av3.extent, [=](concurrency::index<1> idx)  restrict(amp)
		{
			av3[idx] = av1[idx] + av2[idx];
		});
		av3.synchronize();
	}
	catch (const Concurrency::runtime_exception& ex)
	{
		MessageBoxA(NULL, ex.what(), "Error", MB_ICONERROR);
	}
	// Stop timing
	the_amp_clock::time_point end = the_amp_clock::now();
	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Adding the vectors using AMP (data transfer and compute) takes " << time_taken << " ms." << endl;
} // vector_add_amp

  //Attempt to display the mandlebrot
void display_mandlebrot()
{
	HWND window; HDC dc; window = GetActiveWindow(); dc = GetDC(window);

	//This code was too slow for updating the screen so i found a more efficient method
	/*for (int j = 0; j<HEIGHT; j++)
	for (int i = 0; i < WIDTH; i++)
	{
	if (image[j][i] == 0x000000)
	{
	SetPixel(dc, i, j, RGB(0, 0, 0));
	}
	else if (image[j][i] == 0xFFFFFF)
	{
	SetPixel(dc, i, j, RGB(255, 255, 255));
	}
	}*/

	// Creating temp bitmap
	HBITMAP map = CreateBitmap(640, // width
		480, // height
		1, // Color Planes,
		8 * 4, // Size of memory for one pixel in bits (in win32 4 bytes = 4*8 bits)
		(void*)image); // pointer to array
					   // Temp HDC to copy picture
	HDC src = CreateCompatibleDC(dc); // hdc - Device context for window, I've got earlier with GetDC(hWnd) or GetDC(NULL);
	SelectObject(src, map); // Inserting picture into our temp HDC
							// Copy image from temp HDC to window
	BitBlt(dc, // Destination
		0,  // x and
		0,  // y - upper-left corner of place, where we'd like to copy
		640, // width of the region
		480, // height
		src, // source
		0,   // x and
		0,   // y of upper left corner  of part of the source, from where we'd like to copy
		SRCCOPY); // Defined DWORD to juct copy pixels. Watch more on msdn;

	ReleaseDC(window, dc);
	DeleteDC(src);
	DeleteObject(map);
}

void game()
{
	bool finished = false;
	char button = 't';

	float num1 = -2.0;
	float num2 = 1.0;
	float num3 = 1.125;
	float num4 = -1.125;
	float zoom = 1;
	float offsetX = 0;
	float offsetY = 0;

	while (true)
	{
		button = _getch();

		if (button == 'n')
		{
			zoom /= 0.9;
		}

		if (button == 'm')
		{
			zoom *= 0.9;
		}

		if (button == 'w')
		{
			offsetY += 0.5 * zoom;
		}

		if (button == 's')
		{
			offsetY -= 0.5 * zoom;
		}

		if (button == 'a')
		{
			offsetX -= 0.5 * zoom;
		}

		if (button == 'd')
		{
			offsetX += 0.5 * zoom;
		}

		if (button == 'r')
		{
			zoom = 1;
			offsetX = 0;
			offsetY = 0;
		}

		if (button == 'e')
		{
			//exit(0);
			break;
		}

		if (button != NULL)
		{
			compute_mandlebrot((num1 * zoom) + offsetX, (num2 * zoom) + offsetX, (num3 * zoom) + offsetY, (num4 * zoom) + offsetY);
		}

		display_mandlebrot();

		button = NULL;
	}
}

int main(int argc, char *argv[])
{
	// Check AMP support
	//query_AMP_support();

	bool loop = true;
	char input;

	//fill the arrays
	std::vector<int> v1(SIZE, 1);
	std::vector<int> v2(SIZE, 2);
	std::vector<int> v3(SIZE, 0);

	//compare a serial and parallel version of vector addtion
	//time_to_tile(SIZE, v1, v2, v3);
	//vector_add_amp(SIZE, v1, v2, v3);
	//vector_add(SIZE, v1, v2, v3);


	while (loop == true)
	{
		cout << "Welcome to my interactive Mandelbrot Set!" << endl << endl << "Controls:" << endl << "WASD to move the camera, M and N to zoom in and out" << endl << "R to reset back to the default view, E to Exit" << endl << endl << "Press any button to start!" << endl;

		game();

		cout << "Would you like to go back to the Mandlebrot?" << endl;

		cin >> input;

		if (input == 'N' || input == 'n')
		{
			loop = false;
		}
		system("CLS");
	}
	//write_tga("output.tga");

	
	return 0;
} // main