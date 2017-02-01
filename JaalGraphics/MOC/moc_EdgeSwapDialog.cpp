/****************************************************************************
** Meta object code from reading C++ file 'EdgeSwapDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "EdgeSwapDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EdgeSwapDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JEdgeSwapDialog_t {
    QByteArrayData data[9];
    char stringdata0[108];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JEdgeSwapDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JEdgeSwapDialog_t qt_meta_stringdata_JEdgeSwapDialog = {
    {
        QT_MOC_LITERAL(0, 0, 15), // "JEdgeSwapDialog"
        QT_MOC_LITERAL(1, 16, 9), // "startswap"
        QT_MOC_LITERAL(2, 26, 0), // ""
        QT_MOC_LITERAL(3, 27, 6), // "accept"
        QT_MOC_LITERAL(4, 34, 13), // "selectiveEdge"
        QT_MOC_LITERAL(5, 48, 15), // "interactiveEdge"
        QT_MOC_LITERAL(6, 64, 17), // "swapSelectiveEdge"
        QT_MOC_LITERAL(7, 82, 12), // "swapCurrEdge"
        QT_MOC_LITERAL(8, 95, 12) // "skipCurrEdge"

    },
    "JEdgeSwapDialog\0startswap\0\0accept\0"
    "selectiveEdge\0interactiveEdge\0"
    "swapSelectiveEdge\0swapCurrEdge\0"
    "skipCurrEdge"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JEdgeSwapDialog[] = {

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

void JEdgeSwapDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JEdgeSwapDialog *_t = static_cast<JEdgeSwapDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->startswap();
            break;
        case 1:
            _t->accept();
            break;
        case 2:
            _t->selectiveEdge();
            break;
        case 3:
            _t->interactiveEdge();
            break;
        case 4:
            _t->swapSelectiveEdge();
            break;
        case 5:
            _t->swapCurrEdge();
            break;
        case 6:
            _t->skipCurrEdge();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JEdgeSwapDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JEdgeSwapDialog.data,
        qt_meta_data_JEdgeSwapDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JEdgeSwapDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JEdgeSwapDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JEdgeSwapDialog.stringdata0))
        return static_cast<void*>(const_cast< JEdgeSwapDialog*>(this));
    if (!strcmp(_clname, "Ui::EdgeSwapDialog"))
        return static_cast< Ui::EdgeSwapDialog*>(const_cast< JEdgeSwapDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JEdgeSwapDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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
