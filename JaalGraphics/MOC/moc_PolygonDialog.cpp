/****************************************************************************
** Meta object code from reading C++ file 'PolygonDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "PolygonDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'PolygonDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JPolygonDialog_t {
    QByteArrayData data[10];
    char stringdata0[98];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JPolygonDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JPolygonDialog_t qt_meta_stringdata_JPolygonDialog = {
    {
        QT_MOC_LITERAL(0, 0, 14), // "JPolygonDialog"
        QT_MOC_LITERAL(1, 15, 10), // "setPolygon"
        QT_MOC_LITERAL(2, 26, 0), // ""
        QT_MOC_LITERAL(3, 27, 8), // "genShape"
        QT_MOC_LITERAL(4, 36, 11), // "setNumSides"
        QT_MOC_LITERAL(5, 48, 11), // "closeDialog"
        QT_MOC_LITERAL(6, 60, 13), // "keyPressEvent"
        QT_MOC_LITERAL(7, 74, 10), // "QKeyEvent*"
        QT_MOC_LITERAL(8, 85, 1), // "e"
        QT_MOC_LITERAL(9, 87, 10) // "setGenMode"

    },
    "JPolygonDialog\0setPolygon\0\0genShape\0"
    "setNumSides\0closeDialog\0keyPressEvent\0"
    "QKeyEvent*\0e\0setGenMode"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JPolygonDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    6,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    1,       // signalCount

// signals: name, argc, parameters, tag, flags
    1,    0,   44,    2, 0x06 /* Public */,

// slots: name, argc, parameters, tag, flags
    3,    0,   45,    2, 0x08 /* Private */,
    4,    0,   46,    2, 0x08 /* Private */,
    5,    0,   47,    2, 0x08 /* Private */,
    6,    1,   48,    2, 0x08 /* Private */,
    9,    0,   51,    2, 0x08 /* Private */,

// signals: parameters
    QMetaType::Void,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 7,    8,
    QMetaType::Void,

    0        // eod
};

void JPolygonDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JPolygonDialog *_t = static_cast<JPolygonDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->setPolygon();
            break;
        case 1:
            _t->genShape();
            break;
        case 2:
            _t->setNumSides();
            break;
        case 3:
            _t->closeDialog();
            break;
        case 4:
            _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1])));
            break;
        case 5:
            _t->setGenMode();
            break;
        default:
            ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (JPolygonDialog::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&JPolygonDialog::setPolygon)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject JPolygonDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JPolygonDialog.data,
        qt_meta_data_JPolygonDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JPolygonDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JPolygonDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JPolygonDialog.stringdata0))
        return static_cast<void*>(const_cast< JPolygonDialog*>(this));
    if (!strcmp(_clname, "Ui::PolygonDialog"))
        return static_cast< Ui::PolygonDialog*>(const_cast< JPolygonDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JPolygonDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}

// SIGNAL 0
void JPolygonDialog::setPolygon()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
