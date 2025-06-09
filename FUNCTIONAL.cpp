#include <chrono>
#include <cliext/deque>
#include <cliext/utility>
#include <queue>

#include "Windows.h"
#include "framework.h"
#include "pch.h"
enum List { Top, Dawn, Left, Right };
class Interval {
 private:
  std::pair<double, double> start, end;
  double M, R;

 public:
  Interval(std::pair<double, double> _start = std::pair<double, double>(),
           std::pair<double, double> _end = std::pair<double, double>(),
           unsigned short _N = unsigned short()) {
    start = _start, end = _end,
    M = abs(end.second - start.second) /
        pow(abs(end.first - start.first), 1 / double(_N));
  }
  void ChangeCharacteristic(double _m, unsigned short _N) {
    R = pow(end.first - start.first, 1 / double(_N)) +
        (end.second - start.second) * (end.second - start.second) /
            (_m * _m * pow(end.first - start.first, 1 / double(_N))) -
        2 * (end.second + start.second) / _m;
  }
  double GetCharacteristic() { return R; }
  double GetM() { return M; }
  std::pair<double, double> GetEnd() { return end; }
  std::pair<double, double> GetStart() { return start; }
};
class Compare {
 public:
  bool operator()(Interval below, Interval above) {
    if (below.GetCharacteristic() < above.GetCharacteristic()) {
      return true;
    }
    if (below.GetCharacteristic() == above.GetCharacteristic() &&
        below.GetStart().first > above.GetStart().first) {
      return true;
    }
    return false;
  }
};
typedef std::priority_queue<Interval, std::deque<Interval>, Compare>
    Priority_queue;
ref class PeanoCurve_2D {
 private:
  cliext::pair<double, double> ^ x1x2;
  PeanoCurve_2D ^ DawnLeft;
  PeanoCurve_2D ^ DawnRight;
  PeanoCurve_2D ^ TopLeft;
  PeanoCurve_2D ^ TopRight;
  List Type;
  double a, b, c, d;
  unsigned short razvertka;

 public:
  PeanoCurve_2D(unsigned short _razvertka, List _Type, double _a, double _b,
                double _c, double _d) {
    x1x2 = gcnew cliext::pair<double, double>(0.5 * (_a + _b), 0.5 * (_c + _d)),
    Type = _Type, a = _a, b = _b, c = _c, d = _d, razvertka = _razvertka;
    if (_razvertka-- != 0) {
      if (Type == Top) {
        DawnLeft = gcnew PeanoCurve_2D(_razvertka, Right, _a, 0.5 * (_a + _b),
                                       _c, 0.5 * (_c + _d));
        DawnRight = gcnew PeanoCurve_2D(_razvertka, Left, 0.5 * (_a + _b), _b,
                                        _c, 0.5 * (_c + _d));
        TopLeft = gcnew PeanoCurve_2D(_razvertka, Top, _a, 0.5 * (_a + _b),
                                      0.5 * (_c + _d), _d);
        TopRight = gcnew PeanoCurve_2D(_razvertka, Top, 0.5 * (_a + _b), _b,
                                       0.5 * (_c + _d), _d);
      }
      if (Type == Dawn) {
        DawnLeft = gcnew PeanoCurve_2D(_razvertka, Dawn, _a, 0.5 * (_a + _b),
                                       _c, 0.5 * (_c + _d));
        DawnRight = gcnew PeanoCurve_2D(_razvertka, Dawn, 0.5 * (_a + _b), _b,
                                        _c, 0.5 * (_c + _d));
        TopLeft = gcnew PeanoCurve_2D(_razvertka, Right, _a, 0.5 * (_a + _b),
                                      0.5 * (_c + _d), _d);
        TopRight = gcnew PeanoCurve_2D(_razvertka, Left, 0.5 * (_a + _b), _b,
                                       0.5 * (_c + _d), _d);
      }
      if (Type == Left) {
        DawnLeft = gcnew PeanoCurve_2D(_razvertka, Left, _a, 0.5 * (_a + _b),
                                       _c, 0.5 * (_c + _d));
        DawnRight = gcnew PeanoCurve_2D(_razvertka, Top, 0.5 * (_a + _b), _b,
                                        _c, 0.5 * (_c + _d));
        TopLeft = gcnew PeanoCurve_2D(_razvertka, Left, _a, 0.5 * (_a + _b),
                                      0.5 * (_c + _d), _d);
        TopRight = gcnew PeanoCurve_2D(_razvertka, Dawn, 0.5 * (_a + _b), _b,
                                       0.5 * (_c + _d), _d);
      }
      if (Type == Right) {
        DawnLeft = gcnew PeanoCurve_2D(_razvertka, Top, _a, 0.5 * (_a + _b), _c,
                                       0.5 * (_c + _d));
        DawnRight = gcnew PeanoCurve_2D(_razvertka, Right, 0.5 * (_a + _b), _b,
                                        _c, 0.5 * (_c + _d));
        TopLeft = gcnew PeanoCurve_2D(_razvertka, Dawn, _a, 0.5 * (_a + _b),
                                      0.5 * (_c + _d), _d);
        TopRight = gcnew PeanoCurve_2D(_razvertka, Right, 0.5 * (_a + _b), _b,
                                       0.5 * (_c + _d), _d);
      }
    }
  }
  cliext::pair<double, double> HitTest_2D(double x) {
    PeanoCurve_2D ^ tmp = this;
    PeanoCurve_2D ^ Curr = this;
    unsigned short num, i, _razvertka;
    do {
      _razvertka = tmp->razvertka;
      i = 0;
      while (i != _razvertka) {
        num = x * (1 << ++i + i) / (tmp->b - tmp->a);
        if (num == 0) {
          if (Curr->Type == Top || Curr->Type == Right) {
            Curr = Curr->DawnLeft;
          } else {
            Curr = Curr->TopRight;
          }
        }
        if (num == 1) {
          if (Curr->Type == Top || Curr->Type == Left) {
            Curr = Curr->TopLeft;
          } else {
            Curr = Curr->DawnRight;
          }
        }
        if (num == 2) {
          if (Curr->Type == Top || Curr->Type == Right) {
            Curr = Curr->TopRight;
          } else {
            Curr = Curr->DawnLeft;
          }
        }
        if (num == 3) {
          if (Curr->Type == Top || Curr->Type == Left) {
            Curr = Curr->DawnRight;
          } else {
            Curr = Curr->TopLeft;
          }
        }
        x -= num * (tmp->b - tmp->a) * ldexp(1, -i - i);
      }
      tmp = Curr = gcnew PeanoCurve_2D(_razvertka >> 1, Curr->Type, Curr->a,
                                       Curr->b, Curr->c, Curr->d);
    } while (_razvertka != 0);
    return cliext::pair<double, double>(Curr->x1x2);
  }
};
