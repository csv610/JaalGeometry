/****************************************************************************
** Meta object code from reading C++ file 'QuadDiamondsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "QuadDiamondsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QuadDiamondsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JQuadDiamondsDialog_t {
    QByteArrayData data[9];
    char stringdata0[95];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JQuadDiamondsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JQuadDiamondsDialog_t qt_meta_stringdata_JQuadDiamondsDialog = {
    {
        QT_MOC_LITERAL(0, 0, 19), // "JQuadDiamondsDialog"
        QT_MOC_LITERAL(1, 20, 8), // "setColor"
        QT_MOC_LITERAL(2, 29, 0), // ""
        QT_MOC_LITERAL(3, 30, 9), // "removeAll"
        QT_MOC_LITERAL(4, 40, 7), // "meshOpt"
        QT_MOC_LITERAL(5, 48, 17), // "incrementalRemove"
        QT_MOC_LITERAL(6, 66, 14), // "searchDiamonds"
        QT_MOC_LITERAL(7, 81, 6), // "accept"
        QT_MOC_LITERAL(8, 88, 6) // "reject"

    },
    "JQuadDiamondsDialog\0setColor\0\0removeAll\0"
    "meshOpt\0incrementalRemove\0searchDiamonds\0"
    "accept\0reject"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JQuadDiamondsDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    7,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   49,    2, 0x08 /* Private */,
    3,    0,   50,    2, 0x08 /* Private */,
    4,    0,   51,    2, 0x08 /* Private */,
    5,    0,   52,    2, 0x08 /* Private */,
    6,    0,   53,    2, 0x08 /* Private */,
    7,    0,   54,    2, 0x08 /* Private */,
    8,    0,   55,    2, 0x08 /* Private */,

// slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

    0        // eod
};

void JQuadDiamondsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JQuadDiamondsDialog *_t = static_cast<JQuadDiamondsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->setColor();
            break;
        case 1:
            _t->removeAll();
            break;
        case 2:
            _t->meshOpt();
            break;
        case 3:
            _t->incrementalRemove();
            break;
        case 4:
            _t->searchDiamonds();
            break;
        case 5:
            _t->accept();
            break;
        case 6:
            _t->reject();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JQuadDiamondsDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JQuadDiamondsDialog.data,
        qt_meta_data_JQuadDiamondsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JQuadDiamondsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JQuadDiamondsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JQuadDiamondsDialog.stringdata0))
        return static_cast<void*>(const_cast< JQuadDiamondsDialog*>(this));
    if (!strcmp(_clname, "Ui::QuadDiamondsDialog"))
        return static_cast< Ui::QuadDiamondsDialog*>(const_cast< JQuadDiamondsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JQuadDiamondsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
