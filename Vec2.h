//  Real-time Simulation Library
//	Vec2tor.h: Vec2クラス宣言部
//
//  	Copyright(C) 2001,  T. Ishii
//////////////////////////////////////////////////////////////////////
#if !defined(__VEC2_H_INCLUDED)
#define __VEC2_H_INCLUDED 

#include "BaseStream.h"

namespace base {

	//二次元ベクトルクラス
	class Dll Vec2 {
	public:
		Float x, y;
	public:
		//構築
		Vec2();
		//構築
		Vec2(Float _x, Float _y);

		//代入
		void Set(Float _x, Float _y);

		//ゼロにする
		void Zero();

		//要素
		const Float& e(int i) const { return (&x)[i]; }

		//要素
		Float& e(int i) { return (&x)[i]; }

		//要素
		Float& operator[] (int i) { return (&x)[i]; }

		//要素
		const Float& operator[] (int i) const { return (&x)[i]; }

		//四則演算
		Vec2& operator+=(const Vec2 &a);
		Vec2& operator-=(const Vec2 &a);
		Vec2& operator*=(const Float &a);
		Vec2& operator/=(const Float &a);

		//大きさの２乗
		Float SquareMagnitude() const;

		//大きさ
		Float Magnitude() const;

		//正規化
		Vec2 Normalize();

		// 成分の最大を返す
		Float Min() const;

		// 成分の最小を返す
		Float Max() const;

		//絶対値
		Vec2 Abs();

		//ベクトル同士の演算
		friend Vec2 operator+(const Vec2& a, const Vec2& b);
		friend Vec2 operator-(const Vec2& a, const Vec2& b);
		friend Vec2 operator/(const Vec2& a, const Vec2& b);
		// 符合反転
		friend Vec2 operator-(const Vec2& a);

		//スカラーとの演算
		friend Vec2 operator*(const Vec2& a, const Float& b);
		friend Vec2 operator/(const Vec2& a, const Float& b);
		friend Vec2 operator*(const Float& b, const Vec2& a);

		//内積(外積)
		friend Float    Dot(const Vec2& a, const Vec2& b);
		//デターミナント
		friend Float    Det(const Vec2& a, const Vec2& b);
		//外積
		friend Vec2     Cross(const Vec2& a);

		//最大化
		friend Vec2 Minimize(const Vec2& v1, const Vec2& v2);
		//最小化
		friend Vec2 Maximize(const Vec2& v1, const Vec2& v2);

		//比較
		friend int operator == (const Vec2& v1, const Vec2& v2);

		//ストリーム出力
		friend BaseOStream& operator << (BaseOStream& os, const Vec2& a);
		//ストリーム入力
		friend BaseIStream& operator >> (BaseIStream& is, Vec2& a);
	};


	//////////////////////////////////////////////////////////////////////
	//実装部
	//  高速化のためにすべてインライン展開される
	inline
		Vec2::Vec2()
	{
	}

	inline
		Vec2::Vec2(Float _x, Float _y)
	{
		x = _x; y = _y;
	}

	inline void
		Vec2::Set(Float _x, Float _y)
	{
		x = _x; y = _y;
	}

	inline void
		Vec2::Zero()
	{
		x = 0.0; y = 0.0;
	}

	//代入演算
	inline Vec2&
		Vec2::operator+=(const Vec2 &a)
	{
		x = x + a.x; y = y + a.y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator-=(const Vec2 &a)
	{
		x = x - a.x; y = y - a.y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator*=(const Float &a)
	{
		x = a * x; y = a * y;
		return  (*this);
	}

	inline Vec2&
		Vec2::operator/=(const Float &a)
	{
		x = x / a; y = y / a;
		return  (*this);
	}

	//大きさ
	inline Float
		Vec2::SquareMagnitude() const
	{
		return x * x + y * y;
	}

	inline Float
		Vec2::Magnitude()  const
	{
		return sqrt(SquareMagnitude());
	}


	//正規化
	inline Vec2
		Vec2::Normalize()
	{
		Float len = Magnitude();
		Vec2 nor;
		nor.x = x / len;
		nor.y = y / len;
		return nor;
	}

	// 成分の最大,最小を返す
	inline Float
		Vec2::Min()  const
	{
		Float ret = x;
		if (y < ret) ret = y;
		return ret;
	}

	inline Float
		Vec2::Max()  const
	{
		Float ret = x;
		if (ret < y) ret = y;
		return ret;
	}

	//正規化
	inline Vec2
		Vec2::Abs()
	{
		Vec2 r;
		r.x = fabs(x);
		r.y = fabs(y);
		return r;
	}

	///////////////////////////////////////////////////////////////////
	//フレンド関数

	//ベクトル演算
	inline Vec2
		operator+(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret += b;

		return ret;
	}

	inline Vec2
		operator-(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret -= b;

		return ret;
	}

	inline Vec2
		operator/(const Vec2& a, const Vec2& b)
	{
		Vec2 ret = a;
		ret -= b;

		return Vec2(a.x / b.x, a.y / b.y);
	}

	// 符合反転
	inline Vec2
		operator-(const Vec2& a)
	{
		return Vec2(-a.x, -a.y);
	}

	//スカラーとの演算
	inline Vec2
		operator*(const Vec2& a, const Float& b)
	{
		Vec2 ret = a;
		ret *= b;

		return ret;
	}

	inline Vec2
		operator/(const Vec2& a, const Float& b)
	{
		Vec2 ret = a;
		ret /= b;

		return ret;
	}

	inline Vec2
		operator*(const Float& b, const Vec2& a)
	{
		Vec2 ret = a;
		ret *= b;

		return ret;
	}

	//内積
	inline Float
		Dot(const Vec2& a, const Vec2& b)
	{
		return a.x*b.x + a.y*b.y;
	}

	//デターミナント
	inline Float
		Det(const Vec2& a, const Vec2& b)
	{
		return a.x*b.y - a.y*b.x;
	}

	//外積
	inline Vec2
		Cross(const Vec2& a)
	{
		return Vec2(a.y, -a.x);
	}


	//最小化
	inline Vec2
		Minimize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x < v2.x ? v1.x : v2.x,
			v1.y < v2.y ? v1.y : v2.y);
	}

	//最大化
	inline Vec2
		Maximize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x >= v2.x ? v1.x : v2.x,
			v1.y >= v2.y ? v1.y : v2.y);
	}

	//比較
	inline int
		operator == (const Vec2& v1, const Vec2& v2)
	{
		return v1.x == v2.x && v1.y == v2.y;
	}

	//ストリーム
	inline
		BaseOStream& operator << (BaseOStream& os, const Vec2& a)
	{
		if (os.ascii()) {
			os << a.x << " " << a.y;
		}
		else {
			os.write((BYTE*)&a, sizeof(Vec2));
		}
		return os;
	}

	inline
		BaseIStream& operator >> (BaseIStream& is, Vec2& a)
	{
		if (is.ascii()) {
			char str[256];
			is.get_token(str, 256);
			a.x = static_cast<float>(atof(str));
			is.get_token(str, 256);
			a.y = static_cast<float>(atof(str));
		}
		else {
			is.read((BYTE*)&a, sizeof(Vec2));
		}
		return is;
	}

}
#endif/* Not def: __RSVec2TOR_H_INCLUDED */


#pragma once
