//  Real-time Simulation Library
//	Vec2tor.h: Vec2���饹���Բ�
//
//  	Copyright(C) 2001,  T. Ishii
//////////////////////////////////////////////////////////////////////
#if !defined(__VEC2_H_INCLUDED)
#define __VEC2_H_INCLUDED 

#include "BaseStream.h"

namespace base {

	//����Ԫ�٥��ȥ륯�饹
	class Dll Vec2 {
	public:
		Float x, y;
	public:
		//���B
		Vec2();
		//���B
		Vec2(Float _x, Float _y);

		//����
		void Set(Float _x, Float _y);

		//����ˤ���
		void Zero();

		//Ҫ��
		const Float& e(int i) const { return (&x)[i]; }

		//Ҫ��
		Float& e(int i) { return (&x)[i]; }

		//Ҫ��
		Float& operator[] (int i) { return (&x)[i]; }

		//Ҫ��
		const Float& operator[] (int i) const { return (&x)[i]; }

		//�Ąt����
		Vec2& operator+=(const Vec2 &a);
		Vec2& operator-=(const Vec2 &a);
		Vec2& operator*=(const Float &a);
		Vec2& operator/=(const Float &a);

		//�󤭤��Σ��\
		Float SquareMagnitude() const;

		//�󤭤�
		Float Magnitude() const;

		//��Ҏ��
		Vec2 Normalize();

		// �ɷ֤����򷵤�
		Float Min() const;

		// �ɷ֤���С�򷵤�
		Float Max() const;

		//�~����
		Vec2 Abs();

		//�٥��ȥ�ͬʿ������
		friend Vec2 operator+(const Vec2& a, const Vec2& b);
		friend Vec2 operator-(const Vec2& a, const Vec2& b);
		friend Vec2 operator/(const Vec2& a, const Vec2& b);
		// ���Ϸ�ܞ
		friend Vec2 operator-(const Vec2& a);

		//������`�Ȥ�����
		friend Vec2 operator*(const Vec2& a, const Float& b);
		friend Vec2 operator/(const Vec2& a, const Float& b);
		friend Vec2 operator*(const Float& b, const Vec2& a);

		//�ڷe(��e)
		friend Float    Dot(const Vec2& a, const Vec2& b);
		//�ǥ��`�ߥʥ��
		friend Float    Det(const Vec2& a, const Vec2& b);
		//��e
		friend Vec2     Cross(const Vec2& a);

		//���
		friend Vec2 Minimize(const Vec2& v1, const Vec2& v2);
		//��С��
		friend Vec2 Maximize(const Vec2& v1, const Vec2& v2);

		//���^
		friend int operator == (const Vec2& v1, const Vec2& v2);

		//���ȥ�`�����
		friend BaseOStream& operator << (BaseOStream& os, const Vec2& a);
		//���ȥ�`������
		friend BaseIStream& operator >> (BaseIStream& is, Vec2& a);
	};


	//////////////////////////////////////////////////////////////////////
	//�gװ��
	//  ���ٻ��Τ���ˤ��٤ƥ���饤��չ�_�����
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

	//��������
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

	//�󤭤�
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


	//��Ҏ��
	inline Vec2
		Vec2::Normalize()
	{
		Float len = Magnitude();
		Vec2 nor;
		nor.x = x / len;
		nor.y = y / len;
		return nor;
	}

	// �ɷ֤����,��С�򷵤�
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

	//��Ҏ��
	inline Vec2
		Vec2::Abs()
	{
		Vec2 r;
		r.x = fabs(x);
		r.y = fabs(y);
		return r;
	}

	///////////////////////////////////////////////////////////////////
	//�ե����v��

	//�٥��ȥ�����
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

	// ���Ϸ�ܞ
	inline Vec2
		operator-(const Vec2& a)
	{
		return Vec2(-a.x, -a.y);
	}

	//������`�Ȥ�����
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

	//�ڷe
	inline Float
		Dot(const Vec2& a, const Vec2& b)
	{
		return a.x*b.x + a.y*b.y;
	}

	//�ǥ��`�ߥʥ��
	inline Float
		Det(const Vec2& a, const Vec2& b)
	{
		return a.x*b.y - a.y*b.x;
	}

	//��e
	inline Vec2
		Cross(const Vec2& a)
	{
		return Vec2(a.y, -a.x);
	}


	//��С��
	inline Vec2
		Minimize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x < v2.x ? v1.x : v2.x,
			v1.y < v2.y ? v1.y : v2.y);
	}

	//���
	inline Vec2
		Maximize(const Vec2& v1, const Vec2& v2)
	{
		return Vec2(v1.x >= v2.x ? v1.x : v2.x,
			v1.y >= v2.y ? v1.y : v2.y);
	}

	//���^
	inline int
		operator == (const Vec2& v1, const Vec2& v2)
	{
		return v1.x == v2.x && v1.y == v2.y;
	}

	//���ȥ�`��
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
