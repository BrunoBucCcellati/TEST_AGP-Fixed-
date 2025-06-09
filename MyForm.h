#include "FUNCTIONAL.cpp"
namespace TESTAGP {
using namespace System;
using namespace System::Windows::Forms;
public
ref class MyForm : public System::Windows::Forms::Form {
 public:
  double Sign(double Value) {
    if (Value > 0.) {
      return 1;
    }
    return -1;
  }
  double Shag(double _m, double x1, double x2, double y1, double y2,
              unsigned short _N, double _r) {
    if (_N == 1) {
      return 0.5 * (x1 + x2) - 0.5 * (y2 - y1) / _m;
    }
    return 0.5 * (x1 + x2) -
           Sign(y2 - y1) * 0.5 * _r * (y2 - y1) * (y2 - y1) / (_m * _m);
  }
  cliext::deque<double> Base_LNA_1_2_Mer_AGP(
      time_t now, bool mode, unsigned short N, double b, PeanoCurve_2D ^ Curve,
      PeanoCurve_2D ^ Curve_Minus_PI_Na_Dva, unsigned short r, double epsilon,
      unsigned short global_iterations, unsigned short global_local_iterations,
      double a, double c, double d) {
    std::pair<double, double> start, end, start_Minus_PI_Na_Dva,
        end_Minus_PI_Na_Dva, x_Rmax, y_Rmax, x_Rmax_Minus_PI_Na_Dva,
        y_Rmax_Minus_PI_Na_Dva, pred_i_sled_shag,
        pred_i_sled_shag_Minus_PI_Na_Dva, promejutochnaya_tochka,
        promejutochnaya_tochka_Minus_PI_Na_Dva;
    Interval nachalny_otrezok, nachalny_otrezok_Minus_PI_Na_Dva,
        promejutochny_otrezok, promejutochny_otrezok_Minus_PI_Na_Dva, curr,
        curr1, curr_Minus_PI_Na_Dva, curr1_Minus_PI_Na_Dva;
    double Mmax, Mmax_Minus_PI_Na_Dva, m, m_Minus_PI_Na_Dva, dmax,
        dmax_Minus_PI_Na_Dva, eta_shtrih;
    Priority_queue R, R_Minus_PI_Na_Dva, R1, R1_Minus_PI_Na_Dva;
    pred_i_sled_shag = std::pair<double, double>(a, b);
    cliext::deque<double> Extr;
    unsigned short schetchick = 0;
    if (N == 1) {
      HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll");
      typedef double (*sh)(double, time_t);
      sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc");
      start = std::pair<double, double>(a, ShekelFunc(a, now)),
      end = std::pair<double, double>(b, ShekelFunc(b, now)),
      nachalny_otrezok = Interval(start, end, N),
      Mmax = nachalny_otrezok.GetM(), m = r * Mmax,
      x_Rmax = std::pair<double, double>(start.first, end.first),
      y_Rmax = std::pair<double, double>(start.second, end.second),
      R.push(nachalny_otrezok);
      while (true) {
        pred_i_sled_shag.first = pred_i_sled_shag.second,
        promejutochnaya_tochka.first = pred_i_sled_shag.second = Shag(
            m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N, r),
        promejutochnaya_tochka.second =
            ShekelFunc(pred_i_sled_shag.second, now);
        double min = promejutochnaya_tochka.second;
        if (Extr.empty() == false) {
          if (min > Extr.back()) {
            min = Extr.back();
          }
        }
        Extr.push_back(min);
        if (schetchick == global_iterations) {
          return Extr;
        }
        promejutochny_otrezok = R.top(),
        curr = Interval(promejutochny_otrezok.GetStart(),
                        promejutochnaya_tochka, N),
        curr1 =
            Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N),
        R.pop();
        if (mode == true && schetchick > global_local_iterations &&
            schetchick % 2 == 0 && schetchick < 50) {
          if ((std::max)(curr.GetM(), curr1.GetM()) > Mmax ||
              min == promejutochnaya_tochka.second) {
            Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax;
          }
          promejutochny_otrezok = R.top(), R.pop(),
          eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()),
                                  promejutochny_otrezok.GetM()),
          promejutochny_otrezok.ChangeCharacteristic(
              r * Mmax *
                      (promejutochny_otrezok.GetEnd().first -
                       promejutochny_otrezok.GetStart().first) /
                      dmax +
                  eta_shtrih -
                  Mmax *
                      (promejutochny_otrezok.GetEnd().first -
                       promejutochny_otrezok.GetStart().first) /
                      dmax,
              N),
          R.push(promejutochny_otrezok);
          while (R.empty() == false) {
            promejutochny_otrezok = R.top(), R.pop(),
            eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()),
                                    promejutochny_otrezok.GetM()),
            promejutochny_otrezok.ChangeCharacteristic(
                r * Mmax *
                        (promejutochny_otrezok.GetEnd().first -
                         promejutochny_otrezok.GetStart().first) /
                        dmax +
                    eta_shtrih -
                    Mmax *
                        (promejutochny_otrezok.GetEnd().first -
                         promejutochny_otrezok.GetStart().first) /
                        dmax,
                N),
            R1.push(promejutochny_otrezok);
            if (R1.size() == 1) {
              curr.ChangeCharacteristic(
                  r * Mmax * (curr.GetEnd().first - curr.GetStart().first) /
                          dmax +
                      eta_shtrih -
                      Mmax * (curr.GetEnd().first - curr.GetStart().first) /
                          dmax,
                  N),
                  curr1.ChangeCharacteristic(
                      r * Mmax *
                              (curr1.GetEnd().first - curr1.GetStart().first) /
                              dmax +
                          eta_shtrih -
                          Mmax *
                              (curr1.GetEnd().first - curr1.GetStart().first) /
                              dmax,
                      N);
            }
          }
          R = R1, R1 = Priority_queue();
        } else {
          if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax &&
              min != promejutochnaya_tochka.second) {
            curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
          } else {
            Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax,
            curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
            if (mode == true) {
              dmax = (std::max)(
                  pow((curr.GetEnd()).first - (curr.GetStart()).first,
                      (1 / double(N))),
                  pow((curr1.GetEnd()).first - (curr1.GetStart()).first,
                      (1 / double(N))));
            }
            while (R.empty() == false) {
              promejutochny_otrezok = R.top();
              if (mode == true &&
                  pow((promejutochny_otrezok.GetEnd()).first -
                          (promejutochny_otrezok.GetStart()).first,
                      (1 / double(N))) > dmax) {
                dmax = pow((promejutochny_otrezok.GetEnd()).first -
                               (promejutochny_otrezok.GetStart()).first,
                           (1 / double(N)));
              }
              promejutochny_otrezok.ChangeCharacteristic(m, N),
                  R1.push(promejutochny_otrezok), R.pop();
            }
            R = R1, R1 = Priority_queue();
          }
        }
        R.push(curr), R.push(curr1), promejutochny_otrezok = R.top();
        if (abs(promejutochny_otrezok.GetEnd().first -
                promejutochny_otrezok.GetStart().first) < epsilon) {
          return Extr;
        }
        x_Rmax.first = promejutochny_otrezok.GetStart().first,
        x_Rmax.second = promejutochny_otrezok.GetEnd().first,
        y_Rmax.first = promejutochny_otrezok.GetStart().second,
        y_Rmax.second = promejutochny_otrezok.GetEnd().second, schetchick++;
      }
      FreeLibrary(load_function);
    } else {
      HINSTANCE load_function = LoadLibrary(L"TEST_FUNC.dll");
      typedef double (*rr)(double, double);
      rr RastriginFunc = (rr)GetProcAddress(load_function, "RastriginFunc");
      start = std::pair<double, double>(a, RastriginFunc(a, c)),
      end = std::pair<double, double>(b, RastriginFunc(b, c)),
      start_Minus_PI_Na_Dva = std::pair<double, double>(b, RastriginFunc(b, d)),
      end_Minus_PI_Na_Dva = std::pair<double, double>(a, RastriginFunc(a, d)),
      nachalny_otrezok = Interval(start, end, N),
      nachalny_otrezok_Minus_PI_Na_Dva =
          Interval(start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, N),
      Mmax = nachalny_otrezok.GetM(), m = r * Mmax,
      x_Rmax = std::pair<double, double>(start.first, end.first),
      y_Rmax = std::pair<double, double>(start.second, end.second),
      Mmax_Minus_PI_Na_Dva = nachalny_otrezok_Minus_PI_Na_Dva.GetM(),
      m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva, R.push(nachalny_otrezok),
      R_Minus_PI_Na_Dva.push(nachalny_otrezok_Minus_PI_Na_Dva),
      x_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(
          start_Minus_PI_Na_Dva.first, end_Minus_PI_Na_Dva.first),
      y_Rmax_Minus_PI_Na_Dva = std::pair<double, double>(
          start_Minus_PI_Na_Dva.second, end_Minus_PI_Na_Dva.second),
      pred_i_sled_shag_Minus_PI_Na_Dva = std::pair<double, double>(a, b);
      while (true) {
        pred_i_sled_shag.first = pred_i_sled_shag.second,
        promejutochnaya_tochka.first = pred_i_sled_shag.second = Shag(
            m, x_Rmax.first, x_Rmax.second, y_Rmax.first, y_Rmax.second, N, r),
        pred_i_sled_shag_Minus_PI_Na_Dva.first =
            pred_i_sled_shag_Minus_PI_Na_Dva.second,
        promejutochnaya_tochka_Minus_PI_Na_Dva.first =
            pred_i_sled_shag_Minus_PI_Na_Dva.second = Shag(
                m_Minus_PI_Na_Dva, x_Rmax_Minus_PI_Na_Dva.first,
                x_Rmax_Minus_PI_Na_Dva.second, y_Rmax_Minus_PI_Na_Dva.first,
                y_Rmax_Minus_PI_Na_Dva.second, N, r),
        promejutochnaya_tochka.second =
            RastriginFunc(Curve->HitTest_2D(pred_i_sled_shag.second).first,
                          Curve->HitTest_2D(pred_i_sled_shag.second).second),
        promejutochnaya_tochka_Minus_PI_Na_Dva.second = RastriginFunc(
            Curve_Minus_PI_Na_Dva
                ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                .first,
            Curve_Minus_PI_Na_Dva
                ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                .second);
        double min = (std::min)(promejutochnaya_tochka.second,
                                promejutochnaya_tochka_Minus_PI_Na_Dva.second);
        if (Extr.empty() == false) {
          if (min > Extr.back()) {
            min = Extr.back();
          }
        }
        if (min == promejutochnaya_tochka.second) {
          chart2->Series[0]->Points->AddXY(
              Curve->HitTest_2D(pred_i_sled_shag.second).first,
              Curve->HitTest_2D(pred_i_sled_shag.second).second);
        } else {
          chart2->Series[0]->Points->AddXY(
              Curve_Minus_PI_Na_Dva
                  ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                  .first,
              Curve_Minus_PI_Na_Dva
                  ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                  .second);
        }
        Extr.push_back(min);
        if (schetchick == global_iterations) {
          textBox6->Text = Convert::ToString(global_iterations);
          textBox7->Text = Convert::ToString((std::max)(
              abs(promejutochny_otrezok.GetEnd().first -
                  promejutochny_otrezok.GetStart().first),
              abs(promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                  promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first)));
          if (promejutochnaya_tochka.second <
              promejutochnaya_tochka_Minus_PI_Na_Dva.second) {
            textBox4->Text = Convert::ToString(
                Curve->HitTest_2D(pred_i_sled_shag.second).first);
            textBox3->Text = Convert::ToString(
                Curve->HitTest_2D(pred_i_sled_shag.second).second);
          } else {
            textBox4->Text = Convert::ToString(
                Curve_Minus_PI_Na_Dva
                    ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                    .first);
            textBox3->Text = Convert::ToString(
                Curve_Minus_PI_Na_Dva
                    ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                    .second);
          }
          return Extr;
        }
        promejutochny_otrezok = R.top(),
        promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(),
        curr = Interval(promejutochny_otrezok.GetStart(),
                        promejutochnaya_tochka, N),
        curr1 =
            Interval(promejutochnaya_tochka, promejutochny_otrezok.GetEnd(), N),
        R.pop(),
        curr_Minus_PI_Na_Dva =
            Interval(promejutochny_otrezok_Minus_PI_Na_Dva.GetStart(),
                     promejutochnaya_tochka_Minus_PI_Na_Dva, N),
        curr1_Minus_PI_Na_Dva =
            Interval(promejutochnaya_tochka_Minus_PI_Na_Dva,
                     promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd(), N),
        R_Minus_PI_Na_Dva.pop();
        if (mode == true && schetchick > global_local_iterations &&
            schetchick % 2 == 0 && schetchick < 210) {
          if ((std::max)(curr.GetM(), curr1.GetM()) > Mmax ||
              min == promejutochnaya_tochka.second) {
            Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax;
          }
          promejutochny_otrezok = R.top(), R.pop(),
          eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()),
                                  promejutochny_otrezok.GetM()),
          promejutochny_otrezok.ChangeCharacteristic(
              r * Mmax *
                      (promejutochny_otrezok.GetEnd().first -
                       promejutochny_otrezok.GetStart().first) /
                      dmax +
                  eta_shtrih -
                  Mmax *
                      (promejutochny_otrezok.GetEnd().first -
                       promejutochny_otrezok.GetStart().first) /
                      dmax,
              N),
          R.push(promejutochny_otrezok);
          while (R.empty() == false) {
            promejutochny_otrezok = R.top(), R.pop(),
            eta_shtrih = (std::max)((std::max)(curr.GetM(), curr1.GetM()),
                                    promejutochny_otrezok.GetM()),
            promejutochny_otrezok.ChangeCharacteristic(
                r * Mmax *
                        (promejutochny_otrezok.GetEnd().first -
                         promejutochny_otrezok.GetStart().first) /
                        dmax +
                    eta_shtrih -
                    Mmax *
                        (promejutochny_otrezok.GetEnd().first -
                         promejutochny_otrezok.GetStart().first) /
                        dmax,
                N),
            R1.push(promejutochny_otrezok);
            if (R1.size() == 1) {
              curr.ChangeCharacteristic(
                  r * Mmax * (curr.GetEnd().first - curr.GetStart().first) /
                          dmax +
                      eta_shtrih -
                      Mmax * (curr.GetEnd().first - curr.GetStart().first) /
                          dmax,
                  N),
                  curr1.ChangeCharacteristic(
                      r * Mmax *
                              (curr1.GetEnd().first - curr1.GetStart().first) /
                              dmax +
                          eta_shtrih -
                          Mmax *
                              (curr1.GetEnd().first - curr1.GetStart().first) /
                              dmax,
                      N);
            }
          }
          R = R1, R1 = Priority_queue();
          if ((std::max)(curr_Minus_PI_Na_Dva.GetM(),
                         curr1_Minus_PI_Na_Dva.GetM()) > Mmax_Minus_PI_Na_Dva ||
              min == promejutochnaya_tochka_Minus_PI_Na_Dva.second) {
            Mmax_Minus_PI_Na_Dva = (std::max)(curr_Minus_PI_Na_Dva.GetM(),
                                              curr1_Minus_PI_Na_Dva.GetM()),
            m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva;
          }
          promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(),
          R_Minus_PI_Na_Dva.pop(),
          eta_shtrih = (std::max)((std::max)(curr_Minus_PI_Na_Dva.GetM(),
                                             curr1_Minus_PI_Na_Dva.GetM()),
                                  promejutochny_otrezok_Minus_PI_Na_Dva.GetM()),
          promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(
              r * Mmax_Minus_PI_Na_Dva *
                      (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                       promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) /
                      dmax_Minus_PI_Na_Dva +
                  eta_shtrih -
                  Mmax_Minus_PI_Na_Dva *
                      (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                       promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) /
                      dmax_Minus_PI_Na_Dva,
              N),
          R_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva);
          while (R_Minus_PI_Na_Dva.empty() == false) {
            promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(),
            R_Minus_PI_Na_Dva.pop(),
            eta_shtrih =
                (std::max)((std::max)(curr_Minus_PI_Na_Dva.GetM(),
                                      curr1_Minus_PI_Na_Dva.GetM()),
                           promejutochny_otrezok_Minus_PI_Na_Dva.GetM()),
            promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(
                r * Mmax_Minus_PI_Na_Dva *
                        (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                         promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()
                             .first) /
                        dmax_Minus_PI_Na_Dva +
                    eta_shtrih -
                    Mmax_Minus_PI_Na_Dva *
                        (promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                         promejutochny_otrezok_Minus_PI_Na_Dva.GetStart()
                             .first) /
                        dmax_Minus_PI_Na_Dva,
                N),
            R1_Minus_PI_Na_Dva.push(promejutochny_otrezok_Minus_PI_Na_Dva);
            if (R1_Minus_PI_Na_Dva.size() == 1) {
              curr_Minus_PI_Na_Dva.ChangeCharacteristic(
                  r * Mmax_Minus_PI_Na_Dva *
                          (curr_Minus_PI_Na_Dva.GetEnd().first -
                           curr_Minus_PI_Na_Dva.GetStart().first) /
                          dmax_Minus_PI_Na_Dva +
                      eta_shtrih -
                      Mmax_Minus_PI_Na_Dva *
                          (curr_Minus_PI_Na_Dva.GetEnd().first -
                           curr_Minus_PI_Na_Dva.GetStart().first) /
                          dmax_Minus_PI_Na_Dva,
                  N),
                  curr1_Minus_PI_Na_Dva.ChangeCharacteristic(
                      r * Mmax_Minus_PI_Na_Dva *
                              (curr1_Minus_PI_Na_Dva.GetEnd().first -
                               curr1_Minus_PI_Na_Dva.GetStart().first) /
                              dmax_Minus_PI_Na_Dva +
                          eta_shtrih -
                          Mmax_Minus_PI_Na_Dva *
                              (curr1_Minus_PI_Na_Dva.GetEnd().first -
                               curr1_Minus_PI_Na_Dva.GetStart().first) /
                              dmax_Minus_PI_Na_Dva,
                      N);
            }
          }
          R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva,
          R1_Minus_PI_Na_Dva = Priority_queue();
        } else {
          if ((std::max)(curr.GetM(), curr1.GetM()) < Mmax &&
              min != promejutochnaya_tochka.second) {
            curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
          } else {
            Mmax = (std::max)(curr.GetM(), curr1.GetM()), m = r * Mmax,
            curr.ChangeCharacteristic(m, N), curr1.ChangeCharacteristic(m, N);
            if (mode == true) {
              dmax = (std::max)(
                  pow((curr.GetEnd()).first - (curr.GetStart()).first,
                      (1 / double(N))),
                  pow((curr1.GetEnd()).first - (curr1.GetStart()).first,
                      (1 / double(N))));
            }
            while (R.empty() == false) {
              promejutochny_otrezok = R.top(), R.pop();
              if (mode == true &&
                  pow((promejutochny_otrezok.GetEnd()).first -
                          (promejutochny_otrezok.GetStart()).first,
                      (1 / double(N))) > dmax) {
                dmax = pow((promejutochny_otrezok.GetEnd()).first -
                               (promejutochny_otrezok.GetStart()).first,
                           (1 / double(N)));
              }
              promejutochny_otrezok.ChangeCharacteristic(m, N),
                  R1.push(promejutochny_otrezok);
            }
            R = R1, R1 = Priority_queue();
          }
          if ((std::max)(curr_Minus_PI_Na_Dva.GetM(),
                         curr1_Minus_PI_Na_Dva.GetM()) < Mmax_Minus_PI_Na_Dva &&
              min != promejutochnaya_tochka_Minus_PI_Na_Dva.second) {
            curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N),
                curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva,
                                                           N);
          } else {
            Mmax_Minus_PI_Na_Dva = (std::max)(curr_Minus_PI_Na_Dva.GetM(),
                                              curr1_Minus_PI_Na_Dva.GetM()),
            m_Minus_PI_Na_Dva = r * Mmax_Minus_PI_Na_Dva,
            curr_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N),
            curr1_Minus_PI_Na_Dva.ChangeCharacteristic(m_Minus_PI_Na_Dva, N);
            if (mode == true) {
              dmax_Minus_PI_Na_Dva =
                  (std::max)(pow((curr_Minus_PI_Na_Dva.GetEnd()).first -
                                     (curr_Minus_PI_Na_Dva.GetStart()).first,
                                 (1 / double(N))),
                             pow((curr1_Minus_PI_Na_Dva.GetEnd()).first -
                                     (curr1_Minus_PI_Na_Dva.GetStart()).first,
                                 (1 / double(N))));
            }
            while (R_Minus_PI_Na_Dva.empty() == false) {
              promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top(),
              R_Minus_PI_Na_Dva.pop();
              if (mode == true &&
                  pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first -
                          (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart())
                              .first,
                      (1 / double(N))) > dmax_Minus_PI_Na_Dva) {
                dmax_Minus_PI_Na_Dva =
                    pow((promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd()).first -
                            (promejutochny_otrezok_Minus_PI_Na_Dva.GetStart())
                                .first,
                        (1 / double(N)));
              }
              promejutochny_otrezok_Minus_PI_Na_Dva.ChangeCharacteristic(
                  m_Minus_PI_Na_Dva, N),
                  R1_Minus_PI_Na_Dva.push(
                      promejutochny_otrezok_Minus_PI_Na_Dva);
            }
            R_Minus_PI_Na_Dva = R1_Minus_PI_Na_Dva,
            R1_Minus_PI_Na_Dva = Priority_queue();
          }
        }
        R.push(curr), R.push(curr1),
            promejutochny_otrezok = R.top(),
            R_Minus_PI_Na_Dva.push(curr_Minus_PI_Na_Dva),
            R_Minus_PI_Na_Dva.push(curr1_Minus_PI_Na_Dva),
            promejutochny_otrezok_Minus_PI_Na_Dva = R_Minus_PI_Na_Dva.top();
        if (abs(promejutochny_otrezok.GetEnd().first -
                promejutochny_otrezok.GetStart().first) < epsilon &&
            abs(promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first) <
                epsilon) {
          textBox6->Text = Convert::ToString(schetchick);
          textBox7->Text = Convert::ToString((std::max)(
              abs(promejutochny_otrezok.GetEnd().first -
                  promejutochny_otrezok.GetStart().first),
              abs(promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first -
                  promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first)));
          if (promejutochnaya_tochka.second <
              promejutochnaya_tochka_Minus_PI_Na_Dva.second) {
            textBox4->Text = Convert::ToString(
                Curve->HitTest_2D(pred_i_sled_shag.second).first);
            textBox3->Text = Convert::ToString(
                Curve->HitTest_2D(pred_i_sled_shag.second).second);
          } else {
            textBox4->Text = Convert::ToString(
                Curve_Minus_PI_Na_Dva
                    ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                    .first);
            textBox3->Text = Convert::ToString(
                Curve_Minus_PI_Na_Dva
                    ->HitTest_2D(pred_i_sled_shag_Minus_PI_Na_Dva.second)
                    .second);
          }
          return Extr;
        }
        x_Rmax.first = promejutochny_otrezok.GetStart().first,
        x_Rmax.second = promejutochny_otrezok.GetEnd().first,
        y_Rmax.first = promejutochny_otrezok.GetStart().second,
        y_Rmax.second = promejutochny_otrezok.GetEnd().second,
        x_Rmax_Minus_PI_Na_Dva.first =
            promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().first,
        x_Rmax_Minus_PI_Na_Dva.second =
            promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().first,
        y_Rmax_Minus_PI_Na_Dva.first =
            promejutochny_otrezok_Minus_PI_Na_Dva.GetStart().second,
        y_Rmax_Minus_PI_Na_Dva.second =
            promejutochny_otrezok_Minus_PI_Na_Dva.GetEnd().second,
        schetchick++;
      }
      FreeLibrary(load_function);
    }
    return Extr;
  }
  MyForm(void) { InitializeComponent(); }

 protected:
  ~MyForm() {
    if (components) {
      delete components;
    }
  }
  System::Windows::Forms::Button ^ button1;
  System::Windows::Forms::TextBox ^ textBox2;
  System::Windows::Forms::DataVisualization::Charting::Chart ^ chart2;
  System::Windows::Forms::Label ^ label2;
  System::Windows::Forms::TextBox ^ textBox1;
  System::Windows::Forms::TextBox ^ textBox3;
  System::Windows::Forms::TextBox ^ textBox4;
  System::Windows::Forms::TextBox ^ textBox5;
  System::Windows::Forms::TextBox ^ textBox6;
  System::Windows::Forms::Label ^ label6;
  System::Windows::Forms::Label ^ label7;
  System::Windows::Forms::Label ^ label8;
  System::Windows::Forms::Label ^ label9;
  System::Windows::Forms::Label ^ label10;
  System::Windows::Forms::Label ^ label1;
  System::Windows::Forms::TextBox ^ textBox7;
  System::ComponentModel::Container ^ components;
#pragma region Windows Form Designer generated code
  void InitializeComponent(void) {
    System::Windows::Forms::DataVisualization::Charting::ChartArea ^
        chartArea1 =
        (gcnew
             System::Windows::Forms::DataVisualization::Charting::ChartArea());
    System::Windows::Forms::DataVisualization::Charting::Legend ^ legend1 =
        (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
    System::Windows::Forms::DataVisualization::Charting::Series ^ series1 =
        (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
    System::Windows::Forms::DataVisualization::Charting::Series ^ series2 =
        (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
    this->button1 = (gcnew System::Windows::Forms::Button());
    this->textBox2 = (gcnew System::Windows::Forms::TextBox());
    this->chart2 =
        (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
    this->label2 = (gcnew System::Windows::Forms::Label());
    this->textBox1 = (gcnew System::Windows::Forms::TextBox());
    this->textBox3 = (gcnew System::Windows::Forms::TextBox());
    this->textBox4 = (gcnew System::Windows::Forms::TextBox());
    this->textBox5 = (gcnew System::Windows::Forms::TextBox());
    this->textBox6 = (gcnew System::Windows::Forms::TextBox());
    this->label6 = (gcnew System::Windows::Forms::Label());
    this->label7 = (gcnew System::Windows::Forms::Label());
    this->label8 = (gcnew System::Windows::Forms::Label());
    this->label9 = (gcnew System::Windows::Forms::Label());
    this->label10 = (gcnew System::Windows::Forms::Label());
    this->label1 = (gcnew System::Windows::Forms::Label());
    this->textBox7 = (gcnew System::Windows::Forms::TextBox());
    (cli::safe_cast<System::ComponentModel::ISupportInitialize ^>(this->chart2))
        ->BeginInit();
    this->SuspendLayout();
    //
    // button1
    //
    this->button1->BackColor = System::Drawing::SystemColors::Info;
    this->button1->Cursor = System::Windows::Forms::Cursors::Hand;
    this->button1->FlatAppearance->BorderColor =
        System::Drawing::Color::FromArgb(
            static_cast<System::Int32>(static_cast<System::Byte>(64)),
            static_cast<System::Int32>(static_cast<System::Byte>(64)),
            static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->button1->FlatAppearance->BorderSize = 3;
    this->button1->FlatAppearance->MouseDownBackColor =
        System::Drawing::Color::FromArgb(
            static_cast<System::Int32>(static_cast<System::Byte>(128)),
            static_cast<System::Int32>(static_cast<System::Byte>(128)),
            static_cast<System::Int32>(static_cast<System::Byte>(255)));
    this->button1->FlatAppearance->MouseOverBackColor =
        System::Drawing::Color::FromArgb(
            static_cast<System::Int32>(static_cast<System::Byte>(192)),
            static_cast<System::Int32>(static_cast<System::Byte>(192)),
            static_cast<System::Int32>(static_cast<System::Byte>(255)));
    this->button1->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->button1->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 14.25F, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->button1->ForeColor = System::Drawing::SystemColors::ControlDarkDark;
    this->button1->Location = System::Drawing::Point(897, 724);
    this->button1->Name = L"button1";
    this->button1->Size = System::Drawing::Size(275, 75);
    this->button1->TabIndex = 2;
    this->button1->Text = L"SOL";
    this->button1->UseVisualStyleBackColor = false;
    this->button1->Click +=
        gcnew System::EventHandler(this, &MyForm::button1_Click);
    //
    // textBox2
    //
    this->textBox2->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox2->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox2->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox2->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox2->Location = System::Drawing::Point(992, 619);
    this->textBox2->Name = L"textBox2";
    this->textBox2->Size = System::Drawing::Size(180, 29);
    this->textBox2->TabIndex = 4;
    //
    // chart2
    //
    this->chart2->BackColor = System::Drawing::SystemColors::ControlLight;
    chartArea1->AxisX->Interval = 0.1;
    chartArea1->AxisX->IsLabelAutoFit = false;
    chartArea1->AxisX->LabelStyle->Font =
        (gcnew System::Drawing::Font(L"Yu Gothic", 6.75F));
    chartArea1->AxisX->Maximum = 1.8;
    chartArea1->AxisX->Minimum = -2.2;
    chartArea1->AxisX->Title = L"x1";
    chartArea1->AxisX->TitleFont = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 11.25F, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    chartArea1->AxisY->Interval = 0.1;
    chartArea1->AxisY->IsLabelAutoFit = false;
    chartArea1->AxisY->LabelStyle->Font =
        (gcnew System::Drawing::Font(L"Yu Gothic", 7.92F));
    chartArea1->AxisY->Maximum = 1.8;
    chartArea1->AxisY->Minimum = -2.2;
    chartArea1->AxisY->Title = L"x2";
    chartArea1->AxisY->TitleFont = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 11.25F, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    chartArea1->BackColor = System::Drawing::Color::FloralWhite;
    chartArea1->BackGradientStyle = System::Windows::Forms::DataVisualization::
        Charting::GradientStyle::Center;
    chartArea1->BackSecondaryColor = System::Drawing::Color::AliceBlue;
    chartArea1->InnerPlotPosition->Auto = false;
    chartArea1->InnerPlotPosition->Height = 93;
    chartArea1->InnerPlotPosition->Width = 93;
    chartArea1->InnerPlotPosition->X = 3.61F;
    chartArea1->InnerPlotPosition->Y = 1;
    chartArea1->IsSameFontSizeForAllAxes = true;
    chartArea1->Name = L"ChartArea1";
    this->chart2->ChartAreas->Add(chartArea1);
    legend1->BackColor = System::Drawing::Color::Transparent;
    legend1->BackGradientStyle = System::Windows::Forms::DataVisualization::
        Charting::GradientStyle::Center;
    legend1->BackSecondaryColor = System::Drawing::SystemColors::ActiveCaption;
    legend1->BorderColor = System::Drawing::Color::Transparent;
    legend1->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 11.25F, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    legend1->ForeColor = System::Drawing::SystemColors::ActiveCaptionText;
    legend1->HeaderSeparator = System::Windows::Forms::DataVisualization::
        Charting::LegendSeparatorStyle::ThickGradientLine;
    legend1->HeaderSeparatorColor = System::Drawing::Color::IndianRed;
    legend1->IsTextAutoFit = false;
    legend1->ItemColumnSeparator = System::Windows::Forms::DataVisualization::
        Charting::LegendSeparatorStyle::ThickGradientLine;
    legend1->ItemColumnSeparatorColor = System::Drawing::Color::IndianRed;
    legend1->Name = L"Legend1";
    legend1->TableStyle = System::Windows::Forms::DataVisualization::Charting::
        LegendTableStyle::Tall;
    this->chart2->Legends->Add(legend1);
    this->chart2->Location = System::Drawing::Point(12, 14);
    this->chart2->Name = L"chart2";
    series1->BorderWidth = 2;
    series1->ChartArea = L"ChartArea1";
    series1->ChartType = System::Windows::Forms::DataVisualization::Charting::
        SeriesChartType::FastPoint;
    series1->Color = System::Drawing::Color::Blue;
    series1->Legend = L"Legend1";
    series1->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::
        MarkerStyle::Circle;
    series1->Name = L"Точки данных";
    series2->BorderWidth = 2;
    series2->ChartArea = L"ChartArea1";
    series2->ChartType = System::Windows::Forms::DataVisualization::Charting::
        SeriesChartType::FastPoint;
    series2->Color = System::Drawing::Color::Red;
    series2->Legend = L"Legend1";
    series2->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::
        MarkerStyle::Circle;
    series2->Name = L"Точки данных LNA";
    this->chart2->Series->Add(series1);
    this->chart2->Series->Add(series2);
    this->chart2->Size = System::Drawing::Size(1160, 785);
    this->chart2->TabIndex = 5;
    this->chart2->Text = L"chart2";
    //
    // label2
    //
    this->label2->AutoSize = true;
    this->label2->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label2->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label2->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label2->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label2->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label2->ForeColor = System::Drawing::Color::DarkBlue;
    this->label2->Location = System::Drawing::Point(946, 622);
    this->label2->Name = L"label2";
    this->label2->Size = System::Drawing::Size(40, 23);
    this->label2->TabIndex = 8;
    this->label2->Text = L"Extr";
    //
    // textBox1
    //
    this->textBox1->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox1->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox1->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox1->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox1->Location = System::Drawing::Point(992, 654);
    this->textBox1->Name = L"textBox1";
    this->textBox1->Size = System::Drawing::Size(180, 29);
    this->textBox1->TabIndex = 13;
    //
    // textBox3
    //
    this->textBox3->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox3->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox3->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox3->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox3->Location = System::Drawing::Point(992, 584);
    this->textBox3->Name = L"textBox3";
    this->textBox3->Size = System::Drawing::Size(180, 29);
    this->textBox3->TabIndex = 14;
    //
    // textBox4
    //
    this->textBox4->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox4->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox4->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox4->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox4->Location = System::Drawing::Point(992, 549);
    this->textBox4->Name = L"textBox4";
    this->textBox4->Size = System::Drawing::Size(180, 29);
    this->textBox4->TabIndex = 15;
    //
    // textBox5
    //
    this->textBox5->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox5->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox5->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox5->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox5->Location = System::Drawing::Point(992, 514);
    this->textBox5->Name = L"textBox5";
    this->textBox5->Size = System::Drawing::Size(180, 29);
    this->textBox5->TabIndex = 16;
    //
    // textBox6
    //
    this->textBox6->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox6->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox6->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox6->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox6->Location = System::Drawing::Point(992, 479);
    this->textBox6->Name = L"textBox6";
    this->textBox6->Size = System::Drawing::Size(180, 29);
    this->textBox6->TabIndex = 17;
    //
    // label6
    //
    this->label6->AutoSize = true;
    this->label6->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label6->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label6->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label6->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label6->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label6->ForeColor = System::Drawing::Color::DarkBlue;
    this->label6->Location = System::Drawing::Point(911, 657);
    this->label6->Name = L"label6";
    this->label6->Size = System::Drawing::Size(75, 23);
    this->label6->TabIndex = 18;
    this->label6->Text = L"Extr LNA";
    //
    // label7
    //
    this->label7->AutoSize = true;
    this->label7->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label7->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label7->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label7->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label7->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label7->ForeColor = System::Drawing::Color::DarkBlue;
    this->label7->Location = System::Drawing::Point(957, 587);
    this->label7->Name = L"label7";
    this->label7->Size = System::Drawing::Size(29, 23);
    this->label7->TabIndex = 19;
    this->label7->Text = L"x2";
    //
    // label8
    //
    this->label8->AutoSize = true;
    this->label8->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label8->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label8->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label8->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label8->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label8->ForeColor = System::Drawing::Color::DarkBlue;
    this->label8->Location = System::Drawing::Point(960, 552);
    this->label8->Name = L"label8";
    this->label8->Size = System::Drawing::Size(26, 23);
    this->label8->TabIndex = 20;
    this->label8->Text = L"x1";
    //
    // label9
    //
    this->label9->AutoSize = true;
    this->label9->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label9->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label9->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label9->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label9->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label9->ForeColor = System::Drawing::Color::DarkBlue;
    this->label9->Location = System::Drawing::Point(885, 517);
    this->label9->Name = L"label9";
    this->label9->Size = System::Drawing::Size(101, 23);
    this->label9->TabIndex = 21;
    this->label9->Text = L"solving time";
    //
    // label10
    //
    this->label10->AutoSize = true;
    this->label10->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label10->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label10->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label10->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label10->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label10->ForeColor = System::Drawing::Color::DarkBlue;
    this->label10->Location = System::Drawing::Point(903, 482);
    this->label10->Name = L"label10";
    this->label10->Size = System::Drawing::Size(83, 23);
    this->label10->TabIndex = 22;
    this->label10->Text = L"iter count";
    //
    // label1
    //
    this->label1->AutoSize = true;
    this->label1->BackColor = System::Drawing::SystemColors::InactiveCaption;
    this->label1->BorderStyle = System::Windows::Forms::BorderStyle::Fixed3D;
    this->label1->Cursor = System::Windows::Forms::Cursors::Hand;
    this->label1->FlatStyle = System::Windows::Forms::FlatStyle::Flat;
    this->label1->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12,
        static_cast<System::Drawing::FontStyle>(
            ((System::Drawing::FontStyle::Bold |
              System::Drawing::FontStyle::Italic) |
             System::Drawing::FontStyle::Underline)),
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204),
        true));
    this->label1->ForeColor = System::Drawing::Color::DarkBlue;
    this->label1->Location = System::Drawing::Point(911, 692);
    this->label1->Name = L"label1";
    this->label1->Size = System::Drawing::Size(75, 23);
    this->label1->TabIndex = 24;
    this->label1->Text = L"accuracy";
    //
    // textBox7
    //
    this->textBox7->BackColor =
        System::Drawing::SystemColors::ControlLightLight;
    this->textBox7->Cursor = System::Windows::Forms::Cursors::Hand;
    this->textBox7->Font = (gcnew System::Drawing::Font(
        L"Yu Gothic UI", 12, System::Drawing::FontStyle::Bold,
        System::Drawing::GraphicsUnit::Point, static_cast<System::Byte>(204)));
    this->textBox7->ForeColor = System::Drawing::Color::FromArgb(
        static_cast<System::Int32>(static_cast<System::Byte>(64)),
        static_cast<System::Int32>(static_cast<System::Byte>(0)),
        static_cast<System::Int32>(static_cast<System::Byte>(64)));
    this->textBox7->Location = System::Drawing::Point(992, 689);
    this->textBox7->Name = L"textBox7";
    this->textBox7->Size = System::Drawing::Size(180, 29);
    this->textBox7->TabIndex = 23;
    //
    // MyForm
    //
    this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
    this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
    this->ClientSize = System::Drawing::Size(1184, 811);
    this->Controls->Add(this->label1);
    this->Controls->Add(this->textBox7);
    this->Controls->Add(this->label10);
    this->Controls->Add(this->label9);
    this->Controls->Add(this->label8);
    this->Controls->Add(this->label7);
    this->Controls->Add(this->label6);
    this->Controls->Add(this->textBox6);
    this->Controls->Add(this->textBox5);
    this->Controls->Add(this->textBox4);
    this->Controls->Add(this->textBox3);
    this->Controls->Add(this->textBox1);
    this->Controls->Add(this->label2);
    this->Controls->Add(this->textBox2);
    this->Controls->Add(this->button1);
    this->Controls->Add(this->chart2);
    this->Name = L"MyForm";
    this->Text = L"MyForm";
    (cli::safe_cast<System::ComponentModel::ISupportInitialize ^>(this->chart2))
        ->EndInit();
    this->ResumeLayout(false);
    this->PerformLayout();
  }
#pragma endregion
 private:
  PeanoCurve_2D ^ Curve =
      gcnew PeanoCurve_2D(11, List::Top, -2.2, 1.8, -2.2, 1.8);
  PeanoCurve_2D ^ Inverted_Curve =
      gcnew PeanoCurve_2D(11, List::Dawn, -2.2, 1.8, -2.2, 1.8);
  System::Void button1_Click(System::Object ^ sender, System::EventArgs ^ e) {
    chart2->Series[0]->Points->Clear();
    time_t start = std::chrono::duration_cast<std::chrono::nanoseconds>(
                       std::chrono::system_clock::now().time_since_epoch())
                       .count();
    cliext::deque<double> Extr_2D =
        Base_LNA_1_2_Mer_AGP(start, false, 2, 1.8, Curve, Inverted_Curve, 2.5,
                             0.00001, 10000, -1, -2.2, -2.2, 1.8);
    time_t end = std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::system_clock::now().time_since_epoch())
                     .count();
    textBox5->Text =
        Convert::ToString((end - start) / 1'000'000'000.0) + " seconds";
    textBox2->Text = Convert::ToString(Extr_2D.back());
  }
};
}  // namespace TESTAGP
