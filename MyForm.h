#include "FUNCTIONAL.cpp"
namespace TESTAGP {
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		double Sign(double Value)
		{
			if (Value > 0.)
			{
				return 1;
			}
			return -1;
		}
		double Shag(double _m, double x1, double x2, double y1, double y2, unsigned short _N, double _r)
		{
			if (_N == 1)
			{
				return 0.5 * (x1 + x2) - 0.5 * (y2 - y1) / _m;
			}
			return 0.5 * (x1 + x2) - Sign(y2 - y1) * 0.5 * _r * (y2 - y1) * (y2 - y1) / (_m * _m);
		}
		cliext::deque<double> Base_LNA_1_2_Mer_AGP(time_t now, bool mode, unsigned short N, double b, PeanoCurve_2D^ Curve, PeanoCurve_2D^ Curve_Minus_PI_Na_Dva, unsigned short r, double epsilon, unsigned short global_local_iterations, double a, double c, double d)
		{
			std::pair<double, double> start, end, start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, x_Rmax, y_Rmax, x_Rmax_Minus_PI_Na_Dva, y_Rmax_Minus_PI_Na_Dva, pred_i_sled_shag, pred_i_sled_shag_Minus_PI_Na_Dva, promejutochnaya_tochka, promejutochnaya_tochka_Minus_PI_Na_Dva; Interval nachalny_otrezok, nachalny_otrezok_Minus_PI_Na_Dva, promejutochny_otrezok, promejutochny_otrezok_Minus_PI_Na_Dva, curr, curr1, curr_Minus_PI_Na_Dva, curr1_Minus_PI_Na_Dva; double Mmax, Mmax_Minus_PI_Na_Dva, m, m_Minus_PI_Na_Dva, dmax, dmax_Minus_PI_Na_Dva, eta_shtrih; Priority_queue R, R_Minus_PI_Na_Dva, R1, R1_Minus_PI_Na_Dva; pred_i_sled_shag = std::pair<double, double>(a, b); cliext::deque<double> Extr; unsigned short schetchick = 0;
			if (N == 1)
			{
				HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll"); typedef double (*sh) (double, time_t); sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc"); start = std::pair<double, double>(a, ShekelFunc(a, now)), end = std::pair<double, double>(b, ShekelFunc(b, now)), nachalny_otrezok = Interval(start, end, N), Mmax = nachalny_otrezok.GetM(), m = r * Mmax, x_Rmax = std::pair<double, double>(start.first, end.first), y_Rmax = std::pair<double, double>(start.second, end.second), R.push(nachalny_otrezok);
				while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
				{
					pred_i_sled_shag.first = pred_i_sled_shag.second, promejutochnaya_tochka.first = pred_i_sled_shag.second = Shag(m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N, r), promejutochnaya_tochka.second = ShekelFunc(pred_i_sled_shag.second, now); double min = promejutochnaya_tochka.second;
					if (Extr.empty() == false)
					{
						if (min > Extr.back())
						{
							min = Extr.back();
						}
					}
					Extr.push_back(min);
					if (schetchick == 3000)
					{
						return Extr;
					}
					promejutochny_otrezok = R.top(), curr = Interval(promejutochny_otrezok.GetStart(), promejutochnaya_tochka, N), curr1 = Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N), R.pop();
					if (mode == true && schetchick > global_local_iterations && schetchick % 2 == 0 && schetchick < 50)
					{
						if ((std::max)(curr.GetM(), curr1.GetM()) > Mmax || min == promejutochnaya_tochka.second)
						{
							Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax;
						}
						promejutochny_otrezok = R.top(), R.pop(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax, N), R.push(promejutochny_otrezok);
						while (R.empty() == false)
						{
							promejutochny_otrezok = R.top(), R.pop(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax, N), R1.push(promejutochny_otrezok);
							if (R1.size() == 1)
							{
								curr.ChangeCharacteristic(r * Mmax * (curr.GetEnd().first - curr.GetStart().first) / dmax + eta_shtrih - Mmax * (curr.GetEnd().first - curr.GetStart().first) / dmax, N), curr1.ChangeCharacteristic(r * Mmax * (curr1.GetEnd().first - curr1.GetStart().first) / dmax + eta_shtrih - Mmax * (curr1.GetEnd().first - curr1.GetStart().first) / dmax, N);
							}
						}
						R = R1, R1 = Priority_queue();
					}
					else
					{
						if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax && min != promejutochnaya_tochka.second)
						{
							curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
						}
						else
						{
							Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax, curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
							if (mode == true)
							{
								dmax = (std::max)(pow((curr.GetEnd()).first - (curr.GetStart()).first, (1 / double(N))), pow((curr1.GetEnd()).first - (curr1.GetStart()).first, (1 / double(N))));
							}
							while (R.empty() == false)
							{
								promejutochny_otrezok = R.top();
								if (mode == true && pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N))) > dmax)
								{
									dmax = pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N)));
								}
								promejutochny_otrezok.ChangeCharacteristic(m, N), R1.push(promejutochny_otrezok), R.pop();
							}
							R = R1, R1 = Priority_queue();
						}
					}
					R.push(curr), R.push(curr1), promejutochny_otrezok = R.top(), x_Rmax.first = promejutochny_otrezok.GetStart().first, x_Rmax.second = promejutochny_otrezok.GetEnd().first, y_Rmax.first = promejutochny_otrezok.GetStart().second, y_Rmax.second = promejutochny_otrezok.GetEnd().second, schetchick++;
				}
				FreeLibrary(load_function);
			}
			else
			{
				HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll"); typedef double (*grsh) (double, double, time_t); grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc"); start = std::pair<double, double>(a, GrishaginFunc(a, c, now)), end = std::pair<double, double>(b, GrishaginFunc(b, c, now)), start_Minus_PI_Na_Dva = std::pair<double, double>(a, GrishaginFunc(a, c, now)), end_Minus_PI_Na_Dva = std::pair<double, double>(b, GrishaginFunc(a, d, now)), nachalny_otrezok = Interval(start, end, N), nachalny_otrezok_Minus_PI_Na_Dva = Interval(start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, N), Mmax = nachalny_otrezok.GetM(), m = r * Mmax, x_Rmax = std::pair<double, double>(start.first, end.first), y_Rmax = std::pair<double, double>(start.second, end.second), Mmax_Minus_PI_Na_Dva = nachalny_otrezok_Minus_PI_Na_Dva.GetM(), m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva, R.push(nachalny_otrezok), R_Minus_PI_Na_Dva.push(nachalny_otrezok_Minus_PI_Na_Dva), x_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(start_Minus_PI_Na_Dva.first, end_Minus_PI_Na_Dva.first), y_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(start_Minus_PI_Na_Dva.second, end_Minus_PI_Na_Dva.second), pred_i_sled_shag_Minus_PI_Na_Dva = std::pair<double, double>(a, b);
				while (abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon)
				{
					pred_i_sled_shag.first = pred_i_sled_shag.second, promejutochnaya_tochka.first = pred_i_sled_shag.second = Shag(m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N, r), pred_i_sled_shag_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second, promejutochnaya_tochka_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second = Shag(m_Minus_PI_Na_Dva, x_Rmax_Minus_PI_Na_Dva.first, x_Rmax_Minus_PI_Na_Dva.second, y_Rmax_Minus_PI_Na_Dva.first, y_Rmax_Minus_PI_Na_Dva.second, N, r), promejutochnaya_tochka.second = GrishaginFunc(Curve->HitTest_2D(pred_i_sled_shag.second).first, Curve->HitTest_2D(pred_i_sled_shag.second).second, now), promejutochnaya_tochka_Minus_PI_Na_Dva.second = GrishaginFunc(Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).first, Curve_Minus_PI_Na_Dva->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second).second, now); double min = (std::min)(promejutochnaya_tochka.second, promejutochnaya_tochka_Minus_PI_Na_Dva.second);
					if (Extr.empty() == false)
					{
						if (min > Extr.back())
						{
							min = Extr.back();
						}
					}
					Extr.push_back(min);
					if (schetchick == 3000)
					{
						return Extr;
					}
					promejutochny_otrezok = R.top(), promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), curr = Interval(promejutochny_otrezok.GetStart(), promejutochnaya_tochka, N), curr1 = Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N), R.pop(), curr_Minus_PI_Na_Dva = Interval(promejutochny_otrezok_Minus_PI_Na_Dva.GetStart(), promejutochnaya_tochka_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva = Interval(promejutochnaya_tochka_Minus_PI_Na_Dva, promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd(), N), R_Minus_PI_Na_Dva.pop();
					if (mode == true && schetchick > global_local_iterations && schetchick % 2 == 0 && schetchick < 210)
					{
						if ((std::max)(curr.GetM(), curr1.GetM()) > Mmax || min == promejutochnaya_tochka.second)
						{
							Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax;
						}
						promejutochny_otrezok = R.top(), R.pop(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax, N), R.push(promejutochny_otrezok);
						while (R.empty() == false)
						{
							promejutochny_otrezok = R.top(), R.pop(), eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()), promejutochny_otrezok.GetM()), promejutochny_otrezok.ChangeCharacteristic(r * Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax + eta_shtrih - Mmax * (promejutochny_otrezok.GetEnd().first - promejutochny_otrezok.GetStart().first) / dmax, N), R1.push(promejutochny_otrezok);
							if (R1.size() == 1)
							{
								curr.ChangeCharacteristic(r * Mmax * (curr.GetEnd().first - curr.GetStart().first) / dmax + eta_shtrih - Mmax * (curr.GetEnd().first - curr.GetStart().first) / dmax, N), curr1.ChangeCharacteristic(r * Mmax * (curr1.GetEnd().first - curr1.GetStart().first) / dmax + eta_shtrih - Mmax * (curr1.GetEnd().first - curr1.GetStart().first) / dmax, N);
							}
						}
						R = R1, R1 = Priority_queue();
						if ((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()) > Mmax_Minus_PI_Na_Dva || min == promejutochnaya_tochka_Minus_PI_Na_Dva.second)
						{
							Mmax_Minus_PI_Na_Dva = (std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva;
						}
						promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), R_Minus_PI_Na_Dva.pop(), eta_shtrih = (std::max)((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(r* Mmax_Minus_PI_Na_Dva* (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva + eta_shtrih - Mmax_Minus_PI_Na_Dva * (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva, N), R_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva);
						while (R_Minus_PI_Na_Dva.empty() == false)
						{
							promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), R_Minus_PI_Na_Dva.pop(), eta_shtrih = (std::max)((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.GetM()), promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(r * Mmax_Minus_PI_Na_Dva * (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva + eta_shtrih - Mmax_Minus_PI_Na_Dva * (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first - promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva, N), R1_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva);
							if (R1_Minus_PI_Na_Dva.size() == 1)
							{
								curr_Minus_PI_Na_Dva.ChangeCharacteristic(r * Mmax_Minus_PI_Na_Dva * (curr_Minus_PI_Na_Dva.GetEnd().first - curr_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva + eta_shtrih - Mmax_Minus_PI_Na_Dva * (curr_Minus_PI_Na_Dva.GetEnd().first - curr_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva.ChangeCharacteristic(r * Mmax_Minus_PI_Na_Dva * (curr1_Minus_PI_Na_Dva.GetEnd().first - curr1_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva + eta_shtrih - Mmax_Minus_PI_Na_Dva * (curr1_Minus_PI_Na_Dva.GetEnd().first - curr1_Minus_PI_Na_Dva.GetStart().first) / dmax_Minus_PI_Na_Dva, N);
							}
						}
						R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva, R1_Minus_PI_Na_Dva = Priority_queue();
					}
					else
					{
						if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax && min != promejutochnaya_tochka.second)
						{
							curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
						}
						else
						{
							Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax, curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
							if (mode == true)
							{
								dmax = (std::max)(pow((curr.GetEnd()).first - (curr.GetStart()).first, (1 / double(N))), pow((curr1.GetEnd()).first - (curr1.GetStart()).first, (1 / double(N))));
							}
							while (R.empty() == false)
							{
								promejutochny_otrezok = R.top(), R.pop();
								if (mode == true && pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N))) > dmax)
								{
									dmax = pow((promejutochny_otrezok.GetEnd()).first - (promejutochny_otrezok.GetStart()).first, (1 / double(N)));
								}
								promejutochny_otrezok.ChangeCharacteristic(m, N), R1.push(promejutochny_otrezok);
							}
							R = R1, R1 = Priority_queue();
						}
						if ((std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()) < Mmax_Minus_PI_Na_Dva && min != promejutochnaya_tochka_Minus_PI_Na_Dva.second)
						{
							curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N);
						}
						else
						{
							Mmax_Minus_PI_Na_Dva = (std::max)(curr_Minus_PI_Na_Dva.GetM(), curr1_Minus_PI_Na_Dva.GetM()), m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva, curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N);
							if (mode == true)
							{
								dmax_Minus_PI_Na_Dva = (std::max)(pow((curr_Minus_PI_Na_Dva.GetEnd()).first - (curr_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))), pow((curr1_Minus_PI_Na_Dva.GetEnd()).first - (curr1_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))));
							}
							while (R_Minus_PI_Na_Dva.empty() == false)
							{
								promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), R_Minus_PI_Na_Dva.pop();
								if (mode == true && pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first - (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N))) > dmax_Minus_PI_Na_Dva)
								{
									dmax_Minus_PI_Na_Dva = pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first - (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()).first, (1 / double(N)));
								}
								promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N), R1_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva);
							}
							R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva, R1_Minus_PI_Na_Dva = Priority_queue();
						}
					}
					R.push(curr), R.push(curr1), promejutochny_otrezok = R.top(), x_Rmax.first = promejutochny_otrezok.GetStart().first, x_Rmax.second = promejutochny_otrezok.GetEnd().first, y_Rmax.first = promejutochny_otrezok.GetStart().second, y_Rmax.second = promejutochny_otrezok.GetEnd().second, R_Minus_PI_Na_Dva.push(curr_Minus_PI_Na_Dva), R_Minus_PI_Na_Dva.push(curr1_Minus_PI_Na_Dva), promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(), x_Rmax_Minus_PI_Na_Dva.first = promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first, x_Rmax_Minus_PI_Na_Dva.second = promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first, y_Rmax_Minus_PI_Na_Dva.first = promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().second, y_Rmax_Minus_PI_Na_Dva.second = promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().second, schetchick++;
				}
				FreeLibrary(load_function);
			}
			return Extr;
		}
		MyForm(void)
		{
			InitializeComponent();
		}
	protected:
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::DataVisualization::Charting::Chart^ chart1;
	private: System::ComponentModel::Container^ components;
#pragma region Windows Form Designer generated code
		   void InitializeComponent(void)
		   {
			   System::Windows::Forms::DataVisualization::Charting::ChartArea^ chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			   System::Windows::Forms::DataVisualization::Charting::Legend^ legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			   System::Windows::Forms::DataVisualization::Charting::Series^ series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			   System::Windows::Forms::DataVisualization::Charting::Series^ series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			   System::Windows::Forms::DataVisualization::Charting::Series^ series3 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			   System::Windows::Forms::DataVisualization::Charting::Series^ series4 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			   System::Windows::Forms::DataVisualization::Charting::Title^ title1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Title());
			   this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			   (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			   this->SuspendLayout(); 
			   chartArea1->AxisX->Interval = 50;
			   chartArea1->AxisX->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX->IsLabelAutoFit = false;
			   chartArea1->AxisX->LabelAutoFitStyle = static_cast<System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles>(((((((System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::IncreaseFont | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::DecreaseFont)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::StaggeredLabels)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep30)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep45)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep90)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::WordWrap));
			   chartArea1->AxisX->LabelStyle->Font = (gcnew System::Drawing::Font(L"MS Gothic", 9, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   chartArea1->AxisX->LabelStyle->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX->LabelStyle->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX->LineWidth = 3;
			   chartArea1->AxisX->LogarithmBase = 2;
			   chartArea1->AxisX->MajorGrid->LineColor = System::Drawing::Color::Lime;
			   chartArea1->AxisX->Maximum = 3000;
			   chartArea1->AxisX->Minimum = 0;
			   chartArea1->AxisX->MinorGrid->Enabled = true;
			   chartArea1->AxisX->MinorGrid->Interval = 500;
			   chartArea1->AxisX->MinorGrid->LineColor = System::Drawing::Color::Crimson;
			   chartArea1->AxisX->MinorGrid->LineWidth = 4;
			   chartArea1->AxisX->Title = L"NUMBER OF METHOD ITERATIONS";
			   chartArea1->AxisX->TitleFont = (gcnew System::Drawing::Font(L"MS Gothic", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   chartArea1->AxisX2->Interval = 1000;
			   chartArea1->AxisX2->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX2->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisX2->MajorGrid->LineColor = System::Drawing::Color::Red;
			   chartArea1->AxisX2->TitleFont = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8.25F, System::Drawing::FontStyle::Regular,
				   System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
			   chartArea1->AxisY->Interval = 5;
			   chartArea1->AxisY->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisY->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisY->IsLabelAutoFit = false;
			   chartArea1->AxisY->LabelAutoFitStyle = static_cast<System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles>(((((((System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::IncreaseFont | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::DecreaseFont)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::StaggeredLabels)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep30)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep45)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::LabelsAngleStep90)
				   | System::Windows::Forms::DataVisualization::Charting::LabelAutoFitStyles::WordWrap));
			   chartArea1->AxisY->LabelStyle->Font = (gcnew System::Drawing::Font(L"MS Gothic", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   chartArea1->AxisY->LabelStyle->IntervalOffsetType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisY->LabelStyle->IntervalType = System::Windows::Forms::DataVisualization::Charting::DateTimeIntervalType::Number;
			   chartArea1->AxisY->LineWidth = 3;
			   chartArea1->AxisY->MajorGrid->LineColor = System::Drawing::Color::Lime;
			   chartArea1->AxisY->Maximum = 100;
			   chartArea1->AxisY->Minimum = 0;
			   chartArea1->AxisY->MinorGrid->Enabled = true;
			   chartArea1->AxisY->MinorGrid->Interval = 25;
			   chartArea1->AxisY->MinorGrid->LineColor = System::Drawing::Color::Crimson;
			   chartArea1->AxisY->MinorGrid->LineWidth = 4;
			   chartArea1->AxisY->Title = L"PERCENT OF SOLVED";
			   chartArea1->AxisY->TitleFont = (gcnew System::Drawing::Font(L"MS Gothic", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   chartArea1->BackColor = System::Drawing::Color::DimGray;
			   chartArea1->Name = L"ChartArea1";
			   this->chart1->ChartAreas->Add(chartArea1);
			   legend1->BackColor = System::Drawing::Color::SlateGray;
			   legend1->Font = (gcnew System::Drawing::Font(L"MS Gothic", 18, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   legend1->IsTextAutoFit = false;
			   legend1->Name = L"Legend1";
			   this->chart1->Legends->Add(legend1);
			   this->chart1->Location = System::Drawing::Point(12, 12);
			   this->chart1->Name = L"chart1";
			   series1->BorderWidth = 3;
			   series1->ChartArea = L"ChartArea1";
			   series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			   series1->Color = System::Drawing::Color::Red;
			   series1->Font = (gcnew System::Drawing::Font(L"MS Gothic", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   series1->Legend = L"Legend1";
			   series1->MarkerSize = 6;
			   series1->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Square;
			   series1->Name = L"AGP - 2D";
			   series2->BorderWidth = 3;
			   series2->ChartArea = L"ChartArea1";
			   series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			   series2->Color = System::Drawing::Color::Blue;
			   series2->Font = (gcnew System::Drawing::Font(L"MS Gothic", 8.25F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				   static_cast<System::Byte>(204)));
			   series2->Legend = L"Legend1";
			   series2->MarkerSize = 6;
			   series2->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Circle;
			   series2->Name = L"AGP(LNA) - 2D";
			   series3->BorderWidth = 3;
			   series3->ChartArea = L"ChartArea1";
			   series3->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			   series3->Color = System::Drawing::Color::Fuchsia;
			   series3->Legend = L"Legend1";
			   series3->MarkerSize = 6;
			   series3->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Diamond;
			   series3->Name = L"AGP";
			   series4->BorderWidth = 3;
			   series4->ChartArea = L"ChartArea1";
			   series4->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			   series4->Color = System::Drawing::Color::Gold;
			   series4->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 8));
			   series4->Legend = L"Legend1";
			   series4->MarkerSize = 6;
			   series4->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Triangle;
			   series4->Name = L"AGP(LNA)";
			   this->chart1->Series->Add(series1);
			   this->chart1->Series->Add(series2);
			   this->chart1->Series->Add(series3);
			   this->chart1->Series->Add(series4);
			   this->chart1->Size = System::Drawing::Size(2405, 796);
			   this->chart1->TabIndex = 1;
			   this->chart1->Text = L"chart1";
			   title1->BackImageAlignment = System::Windows::Forms::DataVisualization::Charting::ChartImageAlignmentStyle::Center;
			   title1->Font = (gcnew System::Drawing::Font(L"MS Gothic", 18, static_cast<System::Drawing::FontStyle>((System::Drawing::FontStyle::Bold | System::Drawing::FontStyle::Underline)),
				   System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
			   title1->Name = L"Title1";
			   title1->Text = L"OPERATING CHARACTERISTICS OF THE AGP";
			   this->chart1->Titles->Add(title1);
			   this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			   this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			   this->ClientSize = System::Drawing::Size(2429, 820);
			   this->Controls->Add(this->chart1);
			   this->Name = L"MyForm";
			   this->Text = L"MyForm";
			   this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			   (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			   this->ResumeLayout(false);
		   }
#pragma endregion
	private:
		PeanoCurve_2D^ Curve = gcnew PeanoCurve_2D(11, List::Top, 0, 1, 0, 1);
		PeanoCurve_2D^ Curve_Minus_PI_Na_Dva = gcnew PeanoCurve_2D(11, List::Right, 0, 1, 0, 1);
		cli::array<unsigned short>^ procent_correct_2D = gcnew cli::array<unsigned short>(3001);
		cli::array<unsigned short>^ procent_correct_LNA_2D = gcnew cli::array<unsigned short>(3001);
		cli::array<unsigned short>^ procent_correct = gcnew cli::array<unsigned short>(3001);
		cli::array<unsigned short>^ procent_correct_LNA = gcnew cli::array<unsigned short>(3001);
		cli::array<cliext::deque<double>^>^ Extr_2D = gcnew cli::array<cliext::deque<double>^>(1000);
		cli::array<cliext::deque<double>^>^ Extr_LNA_2D = gcnew cli::array<cliext::deque<double>^>(1000);
		cli::array<cliext::deque<double>^>^ Extr = gcnew cli::array<cliext::deque<double>^>(1000);
		cli::array<cliext::deque<double>^>^ Extr_LNA = gcnew cli::array<cliext::deque<double>^>(1000);
		System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e)
		{
			unsigned short num = 0, MIN_LNA_SIZE = 3000, MIN_LNA_SIZE_2D = 3000;
			while (num < 1000)
			{
				time_t now = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
				HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll"); typedef double (*grsh) (double, double, time_t); grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc"); typedef double (*sh) (double, time_t); sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc");
				unsigned short i = 0, j; double MIN_GRSH = 0, MIN_SH = 0;
				while (i < 256)
				{
					j = 0;
					while (j < 256)
					{
						if (GrishaginFunc(ldexp(1, -9) + ldexp(1, -8) * i, ldexp(1, -9) + ldexp(1, -8) * (j += 1), now) < MIN_GRSH)
						{
							MIN_GRSH = GrishaginFunc(ldexp(1, -9) + ldexp(1, -8) * i, ldexp(1, -9) + ldexp(1, -8) * j, now);
						}
					}
					i++;
				}
				i = 2000;
				while (i > 0)
				{
					if (ShekelFunc((i -= 1) * 0.005, now) < MIN_SH)
					{
						MIN_SH = ShekelFunc(i * 0.005, now);
					}
				}
				FreeLibrary(load_function);
				Extr_2D[num] = gcnew cliext::deque<double>(Base_LNA_1_2_Mer_AGP(now, false, 2, 1, Curve, Curve_Minus_PI_Na_Dva, 5, 5 * pow(10, -15), -1, 0, 0, 1));
				Extr_LNA_2D[num] = gcnew cliext::deque<double>(Base_LNA_1_2_Mer_AGP(now, true, 2, 1, Curve, Curve_Minus_PI_Na_Dva, 3, 5 * pow(10, -15), 198, 0, 0, 1));
				Extr[num] = gcnew cliext::deque<double>(Base_LNA_1_2_Mer_AGP(now, false, 1, 10, nullptr, nullptr, 3, 5 * pow(10, -15), -1, 0, 0, 1));
				Extr_LNA[num] = gcnew cliext::deque<double>(Base_LNA_1_2_Mer_AGP(now, true, 1, 10, nullptr, nullptr, 3, 5 * pow(10, -15), 38, 0, 0, 1));
				if (Extr_LNA[num]->size() < MIN_LNA_SIZE)
				{
					MIN_LNA_SIZE = Extr_LNA[num]->size();
				}
				if (Extr_LNA_2D[num]->size() < MIN_LNA_SIZE_2D)
				{
					MIN_LNA_SIZE_2D = Extr_LNA_2D[num]->size();
				}
				while (Extr_2D[num]->empty() == false)
				{
					if (Extr_2D[num]->front() < MIN_GRSH || Extr_2D[num]->front() - MIN_GRSH < 0.01)
					{
						procent_correct_2D[i]++;
					}
					if (num == 999)
					{
						chart1->Series[0]->Points->AddXY(i, procent_correct_2D[i] * 0.1);
					}
					Extr_2D[num]->pop_front();
					if (Extr[num]->empty() == false)
					{
						if (Extr[num]->front() < MIN_SH || Extr[num]->front() - MIN_SH < 0.001)
						{
							procent_correct[i]++;
						}
						if (num == 999)
						{
							chart1->Series[2]->Points->AddXY(i, procent_correct[i] * 0.1);
						}
						Extr[num]->pop_front();
					}
					if (Extr_LNA_2D[num]->empty() == false)
					{
						if (Extr_LNA_2D[num]->front() < MIN_GRSH || Extr_LNA_2D[num]->front() - MIN_GRSH < 0.01)
						{
							procent_correct_LNA_2D[i]++;
						}
						if (num == 999)
						{
							if (i < MIN_LNA_SIZE_2D)
							{
								chart1->Series[1]->Points->AddXY(i, procent_correct_LNA_2D[i] * 0.1);
							}
						}
						Extr_LNA_2D[num]->pop_front();
					}
					if (Extr_LNA[num]->empty() == false)
					{
						if (Extr_LNA[num]->front() < MIN_SH || Extr_LNA[num]->front() - MIN_SH < 0.001)
						{
							procent_correct_LNA[i]++;
						}
						if (num == 999)
						{
							if (i < MIN_LNA_SIZE)
							{
								chart1->Series[3]->Points->AddXY(i, procent_correct_LNA[i] * 0.1);
							}
						}
						Extr_LNA[num]->pop_front();
					}
					i++;
				}
				num++;
			}
		}
};
}