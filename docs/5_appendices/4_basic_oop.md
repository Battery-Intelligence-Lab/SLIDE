---
layout: default
title: OOP basics
nav_order: 4
---


# The basics of object-oriented programming

<p style='text-align: justify;'>
The code is written in an object-oriented programming style. Below is a very short introduction to the basics of object-oriented programming. It should be enough to understand the code from this project, but you can find much more information on the internet if things are not clear.
</p>

## Classes and objects

<p style='text-align: justify;'>
The basic elements are classes and objects. They are best explained with a real-world example. Think of a car. A car has certain properties (e.g. a colour, a brand, a year in which the car was built, etc.), a given state (e.g. the speed, location, etc.). A class is an abstract representation of what ‘a car’ is (i.e. a ‘template’ for a class). The properties and states are represented as so-called ‘class variables’. They are defined somewhere in the class (without giving values for those variables, the class is just the template defining what properties a car has). In pseudo-code, our Car-class looks like:
</p>

```cpp
Class Car {
  string brand;
  string colour;
  int buildYear;
  double speed;
  //...
}
```

Up to here, a class is very similar to a struct in MATLAB, C or other programming languages. The main difference is that in the real-world, there are clearly defined ‘things’ a car can ‘do’ (e.g. accelerate, turn, etc.). While in MATLAB or other non-object-oriented programming functions to ‘do things’ have to be defined externally, they are included in a class in object-oriented code. So a class groups both the variables AND the functions related to the ‘thing’ the class represents. Let’s add a function to accelerate the car to our class. The function has one variable, the increase in speed ds. Our Car-class could for instance look like:

```cpp
Class Car {
  string brand;
  string colour;
  int buildYear;
  double speed;
  //...
  accelerate(double ds) { Speed = speed + ds; }
}
```

Note that you can access the class variables in functions as if they were local variables. The class variables are shared between all functions of the class, so if one function changes the value than the value will also be different in the other functions. In fact, they are some sort of ‘local global’ variables for all the functions inside our Car-class.

An Object is an ‘instance’ of our class, i.e. with concrete values for the class variables. E.g. we can make an object called c1, and set its variables as follows:

```cpp
Car c1 = new Car;
c1.brand = “Opel”;
c1.colour = “green”;
c1.buildYear = 2015;
c1.speed = 0;
```

In this case, we have the Object c1 which is a green Opel car made in 2015 which is stationary. If we now want to accelerate our car by 1 m/s, we can call the function from the object:

```cpp
c1.accelerate(1);
```

Now our green Opel car will be driving at 1 m/s, i.e. `c1.speed` will be 1.

One of the advantages is that it significantly simplifies the book keeping if we have many cars in our simulation. E.g. if we want to make two more cars, c2 and c3, we type the following:

```cpp
Car c2 = new Car;
c2.brand = “Ford”;
c2.colour = “red”;
c2.buildYear = 2016;
c2.speed = 0;

Car c3 = new Car;
c3.brand = “Renault”;
c3.colour = “yellow”;
c3.buildYear = 2010;
c3.speed = 0;
```

<p style='text-align: justify;'>
Now we can use our three cars independently of each other and we don’t have to remember the state and parameters of each car. E.g. if we type c2.accelerate(2), the Ford is driving at 2m/s, while the Opel is still driving at 1m/s and the Renault is still stationary. This means we don’t have to remember which speed corresponded to which car (because each car remembers its own speed).
</p>

## Encapsulation

<p style='text-align: justify;'>
Encapsulation means that inside a class, you ‘hide complexity’ to the outside world. A class has ‘private’ features which are only accessible from within the class and ‘public’ features which can be accessed externally. This reduces complexity and decreases the potential for errors in the code. It is good coding practice to make all class variable private, and make public functions to get their values (‘getters’) and to set their values (‘setters’). In the setters, you can then check if the values are realistic, which will avoid getting errors in the code. The user can only change the value of the variable with the setter, so every time we set a new value, it will be checked on validity, thus avoiding getting wrong values. E.g. if we only allow build years after 1800, our Car could look like:
</p>

```cpp
Class Car {
private:
  string brand;
  string colour;
  int buildYear;
  double speed;

public:
  getBrand() { return brand; }
  getColour() { return colour; }
  getbuildYear() { return buildYear; }
  getSpeed() { return speed; }
  setBrand(string newBrand) { brand = newBrand; }
  setColour(string newColour) { colour = newColour; }
  setbuildYear(int newYear) {
    if
      newYear > 1800 { buildYear = newYear; }
    else {
      throw error;
    }
  }
  setSpeed(int newSpeed) { speed = newSpeed; }
  accelerate(double ds) { speed = speed + ds; }
}
```

## Relation between classes

Another great feature of object-oriented classes is that they can be linked to each other. Suppose we have a class to represent a tire, e.g. with a radius and the air pressure in the tire:

```cpp
Class wheel {
private:
  double radius;
  double pressure;

public:
  getRadius { return radius; }
  getPressure { return pressure; }
  setRadius(int newRadius) { radius = newRadius; }
  setPressure(int newPressure) { pressure = newPressure; }
}
```

Now we can *add the wheels* to our car. Class variables can be of any type, allowing this sort of relations.

```cpp
Class Car {
private:
  string brand;
  string colour;
  int buildYear;
  double speed;
  Wheel leftFront;
  Wheel leftRear;
  Wheel rightFront;
  Wheel rightRear;

public:
  getBrand() { return brand; }
  getColour() { return colour; }
  getbuildYear() { return buildYear; }
  getSpeed() { return speed; }
  getWheels() { return {leftFront, leftRear, rightFront, rightRear}; }
  setBrand(string newBrand) { brand = newBrand; }
  setColour(string newColour) { colour = newColour; }
  setbuildYear(int newYear) {
    if
      newYear > 1800 { buildYear = newYear; }
    else {
      throw error;
    }
  }
  setSpeed(int newSpeed) { speed = newSpeed; }
  setWheels(Wheel newlf, Wheel newlr, Wheel newrf, Wheel newrr) {
    leftFront = newlf;
    leftRear = newlr;
    rightFront = newrf;
    rightRear = newrr;
  }
  accelerate(double ds) { speed = speed + ds; }
}
```

And let us make 4 wheels with the same radius but different pressures:

```cpp
Wheel w1 = new Wheel;
w1.setRadius(0.5);
w1.setPressure(10);
Wheel w2 = new Wheel;
w2.setRadius(0.5);
w2.setPressure(11);
Wheel w3 = new Wheel;
w3.setRadius(0.5);
w3.setPressure(9);
Wheel w3 = new Wheel;
w3.setRadius(0.5);
w3.setPressure(12);
Wheel w4 = new Wheel;
w4.setRadius(0.5);
w4.setPressure(10);
```

Now we can add the wheels to our cars:

```cpp
c1.setWheels(w1, w2, w3, w4);
```

<p style='text-align: justify;'>
We could do the same with other wheels and other cars. As you can see, this will reduce the amount of bookkeeping we have to do to remember which wheel was part of which car (because the cars remember themselves).
Also keeping consistency between our objects becomes easier. Suppose that wheels also have a class-variable rotationSpeed, which gives the speed at which wheels are rotating (and we have implemented the standard getters and setters for the rotation speed in our Wheel-class). Of course, the speed of the car is related to the rotation speed of the wheels (the speed of the wheels is the speed of the car divided by the radius of the wheels). Every time we update the speed of the car, we have to update the speed of the wheels as well. The code below shows how the functions setSpeed would look like when we take advantage of these features. In our Car-class we can access the radius of the wheels using the getter from our Wheel class (getRadius), and set the rotation speed of our wheels using the setter from Wheel (setRotationSpeed(double newRotationSpeed) ).
</p>

```cpp
setSpeed(int newSpeed) {
  speed = newSpeed;
  leftFront.setRotationSpeed(newspeed / leftFront.getRadius());
  leftRear.setRotationSpeed(newspeed / leftRear.getRadius());
  rightFront.setRotationSpeed(newspeed / rightFront.getRadius());
  rightRear.setRotationSpeed(newspeed / rightRear.getRadius());
}
```

Now in the function accelerate, we can simply call the setSpeed function such that we don’t have to worry about the rotation speed of the wheels when accelerating the car:

```cpp
accelerate(double ds) { setSpeed(getSpeed() + ds); }
```

Our Car-class has fully ‘hidden’ the wheels to the outside world. When you change the speed of the car, you don’t have to think about the wheels. This is a good example of the ‘encapsulation’ principle mentioned before.

## And much more

Object-oriented programming has many other features, such as inheritance (if a child class ‘inherits’ from a parent class, it is as if the child class extends the parent class, i.e. it has all the features from the parent class plus any additional featured implemented in the child class).

There are again loads of tutorials online. A short introduction is given on

[http://ee402.eeng.dcu.ie/introduction/chapter-1---introduction-to-object-oriented-programming](http://ee402.eeng.dcu.ie/introduction/chapter-1---introduction-to-object-oriented-programming)