/****************************************************************************
** Meta object code from reading C++ file 'CurveGenDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/CurveGenDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'CurveGenDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JCurveGenDialog_t {
    QByteArrayData data[12];
    char stringdata0[147];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JCurveGenDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JCurveGenDialog_t qt_meta_stringdata_JCurveGenDialog = {
    {
QT_MOC_LITERAL(0, 0, 15), // "JCurveGenDialog"
QT_MOC_LITERAL(1, 16, 10), // "closeCurve"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 11), // "deleteCurve"
QT_MOC_LITERAL(4, 40, 11), // "deletePoint"
QT_MOC_LITERAL(5, 52, 15), // "openKnotsDialog"
QT_MOC_LITERAL(6, 68, 8), // "getKnots"
QT_MOC_LITERAL(7, 77, 14), // "refineSegments"
QT_MOC_LITERAL(8, 92, 14), // "startSketching"
QT_MOC_LITERAL(9, 107, 12), // "createCanvas"
QT_MOC_LITERAL(10, 120, 14), // "setCanvasColor"
QT_MOC_LITERAL(11, 135, 11) // "closeDialog"

    },
    "JCurveGenDialog\0closeCurve\0\0deleteCurve\0"
    "deletePoint\0openKnotsDialog\0getKnots\0"
    "refineSegments\0startSketching\0"
    "createCanvas\0setCanvasColor\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JCurveGenDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   64,    2, 0x08 /* Private */,
       3,    0,   65,    2, 0x08 /* Private */,
       4,    0,   66,    2, 0x08 /* Private */,
       5,    0,   67,    2, 0x08 /* Private */,
       6,    0,   68,    2, 0x08 /* Private */,
       7,    0,   69,    2, 0x08 /* Private */,
       8,    0,   70,    2, 0x08 /* Private */,
       9,    0,   71,    2, 0x08 /* Private */,
      10,    0,   72,    2, 0x08 /* Private */,
      11,    0,   73,    2, 0x08 /* Private */,

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

       0        // eod
};

void JCurveGenDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JCurveGenDialog *_t = static_cast<JCurveGenDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->closeCurve(); break;
        case 1: _t->deleteCurve(); break;
        case 2: _t->deletePoint(); break;
        case 3: _t->openKnotsDialog(); break;
        case 4: _t->getKnots(); break;
        case 5: _t->refineSegments(); break;
        case 6: _t->startSketching(); break;
        case 7: _t->createCanvas(); break;
        case 8: _t->setCanvasColor(); break;
        case 9: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JCurveGenDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JCurveGenDialog.data,
      qt_meta_data_JCurveGenDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JCurveGenDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JCurveGenDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JCurveGenDialog.stringdata0))
        return static_cast<void*>(const_cast< JCurveGenDialog*>(this));
    if (!strcmp(_clname, "Ui::CurveGenDialog"))
        return static_cast< Ui::CurveGenDialog*>(const_cast< JCurveGenDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JCurveGenDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 10)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 10;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 10)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 10;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
