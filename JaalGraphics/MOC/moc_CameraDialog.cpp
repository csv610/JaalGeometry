/****************************************************************************
** Meta object code from reading C++ file 'CameraDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "CameraDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'CameraDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JCameraDialog_t {
    QByteArrayData data[13];
    char stringdata0[176];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JCameraDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JCameraDialog_t qt_meta_stringdata_JCameraDialog = {
    {
QT_MOC_LITERAL(0, 0, 13), // "JCameraDialog"
QT_MOC_LITERAL(1, 14, 11), // "setViewSide"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 12), // "setViewPoint"
QT_MOC_LITERAL(4, 40, 9), // "setCamera"
QT_MOC_LITERAL(5, 50, 13), // "setProjection"
QT_MOC_LITERAL(6, 64, 11), // "entireScene"
QT_MOC_LITERAL(7, 76, 14), // "setSceneRadius"
QT_MOC_LITERAL(8, 91, 21), // "setRotationConstraint"
QT_MOC_LITERAL(9, 113, 24), // "setTranslationConstraint"
QT_MOC_LITERAL(10, 138, 11), // "setPosition"
QT_MOC_LITERAL(11, 150, 13), // "setWindowSize"
QT_MOC_LITERAL(12, 164, 11) // "closeDialog"

    },
    "JCameraDialog\0setViewSide\0\0setViewPoint\0"
    "setCamera\0setProjection\0entireScene\0"
    "setSceneRadius\0setRotationConstraint\0"
    "setTranslationConstraint\0setPosition\0"
    "setWindowSize\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JCameraDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   69,    2, 0x08 /* Private */,
       3,    0,   70,    2, 0x08 /* Private */,
       4,    0,   71,    2, 0x08 /* Private */,
       5,    0,   72,    2, 0x08 /* Private */,
       6,    0,   73,    2, 0x08 /* Private */,
       7,    0,   74,    2, 0x08 /* Private */,
       8,    0,   75,    2, 0x08 /* Private */,
       9,    0,   76,    2, 0x08 /* Private */,
      10,    0,   77,    2, 0x08 /* Private */,
      11,    0,   78,    2, 0x08 /* Private */,
      12,    0,   79,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JCameraDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JCameraDialog *_t = static_cast<JCameraDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setViewSide(); break;
        case 1: _t->setViewPoint(); break;
        case 2: _t->setCamera(); break;
        case 3: _t->setProjection(); break;
        case 4: _t->entireScene(); break;
        case 5: _t->setSceneRadius(); break;
        case 6: _t->setRotationConstraint(); break;
        case 7: _t->setTranslationConstraint(); break;
        case 8: _t->setPosition(); break;
        case 9: _t->setWindowSize(); break;
        case 10: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JCameraDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JCameraDialog.data,
      qt_meta_data_JCameraDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JCameraDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JCameraDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JCameraDialog.stringdata0))
        return static_cast<void*>(const_cast< JCameraDialog*>(this));
    if (!strcmp(_clname, "Ui::CameraDialog"))
        return static_cast< Ui::CameraDialog*>(const_cast< JCameraDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JCameraDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
