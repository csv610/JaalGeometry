/****************************************************************************
** Meta object code from reading C++ file 'GenSimpleShapeDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GenSimpleShapeDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GenSimpleShapeDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JGenSimpleShapeDialog_t {
    QByteArrayData data[13];
    char stringdata0[181];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JGenSimpleShapeDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JGenSimpleShapeDialog_t qt_meta_stringdata_JGenSimpleShapeDialog = {
    {
        QT_MOC_LITERAL(0, 0, 21), // "JGenSimpleShapeDialog"
        QT_MOC_LITERAL(1, 22, 4), // "init"
        QT_MOC_LITERAL(2, 27, 0), // ""
        QT_MOC_LITERAL(3, 28, 8), // "genShape"
        QT_MOC_LITERAL(4, 37, 10), // "getPolygon"
        QT_MOC_LITERAL(5, 48, 15), // "openKnotsDialog"
        QT_MOC_LITERAL(6, 64, 15), // "openCurveDialog"
        QT_MOC_LITERAL(7, 80, 17), // "openLSystemDialog"
        QT_MOC_LITERAL(8, 98, 17), // "openPolygonDialog"
        QT_MOC_LITERAL(9, 116, 19), // "openPolyhedraDialog"
        QT_MOC_LITERAL(10, 136, 20), // "openImageEdgesDialog"
        QT_MOC_LITERAL(11, 157, 11), // "acceptShape"
        QT_MOC_LITERAL(12, 169, 11) // "cancelShape"

    },
    "JGenSimpleShapeDialog\0init\0\0genShape\0"
    "getPolygon\0openKnotsDialog\0openCurveDialog\0"
    "openLSystemDialog\0openPolygonDialog\0"
    "openPolyhedraDialog\0openImageEdgesDialog\0"
    "acceptShape\0cancelShape"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JGenSimpleShapeDialog[] = {

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

void JGenSimpleShapeDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JGenSimpleShapeDialog *_t = static_cast<JGenSimpleShapeDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->init();
            break;
        case 1:
            _t->genShape();
            break;
        case 2:
            _t->getPolygon();
            break;
        case 3:
            _t->openKnotsDialog();
            break;
        case 4:
            _t->openCurveDialog();
            break;
        case 5:
            _t->openLSystemDialog();
            break;
        case 6:
            _t->openPolygonDialog();
            break;
        case 7:
            _t->openPolyhedraDialog();
            break;
        case 8:
            _t->openImageEdgesDialog();
            break;
        case 9:
            _t->acceptShape();
            break;
        case 10:
            _t->cancelShape();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JGenSimpleShapeDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JGenSimpleShapeDialog.data,
        qt_meta_data_JGenSimpleShapeDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JGenSimpleShapeDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JGenSimpleShapeDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JGenSimpleShapeDialog.stringdata0))
        return static_cast<void*>(const_cast< JGenSimpleShapeDialog*>(this));
    if (!strcmp(_clname, "Ui::GenSimpleShapeDialog"))
        return static_cast< Ui::GenSimpleShapeDialog*>(const_cast< JGenSimpleShapeDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JGenSimpleShapeDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
