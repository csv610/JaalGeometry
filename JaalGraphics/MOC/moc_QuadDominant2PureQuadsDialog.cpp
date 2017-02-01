/****************************************************************************
** Meta object code from reading C++ file 'QuadDominant2PureQuadsDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/QuadDominant2PureQuadsDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'QuadDominant2PureQuadsDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JQuadDominant2PureQuadsDialog_t {
    QByteArrayData data[10];
    char stringdata0[125];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JQuadDominant2PureQuadsDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JQuadDominant2PureQuadsDialog_t qt_meta_stringdata_JQuadDominant2PureQuadsDialog = {
    {
QT_MOC_LITERAL(0, 0, 29), // "JQuadDominant2PureQuadsDialog"
QT_MOC_LITERAL(1, 30, 9), // "enumFaces"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 11), // "getNewStrip"
QT_MOC_LITERAL(4, 53, 11), // "remeshStrip"
QT_MOC_LITERAL(5, 65, 12), // "getAllStrips"
QT_MOC_LITERAL(6, 78, 15), // "remeshAllStrips"
QT_MOC_LITERAL(7, 94, 9), // "refineAll"
QT_MOC_LITERAL(8, 104, 8), // "clearAll"
QT_MOC_LITERAL(9, 113, 11) // "closeDialog"

    },
    "JQuadDominant2PureQuadsDialog\0enumFaces\0"
    "\0getNewStrip\0remeshStrip\0getAllStrips\0"
    "remeshAllStrips\0refineAll\0clearAll\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JQuadDominant2PureQuadsDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   54,    2, 0x08 /* Private */,
       3,    0,   55,    2, 0x08 /* Private */,
       4,    0,   56,    2, 0x08 /* Private */,
       5,    0,   57,    2, 0x08 /* Private */,
       6,    0,   58,    2, 0x08 /* Private */,
       7,    0,   59,    2, 0x08 /* Private */,
       8,    0,   60,    2, 0x08 /* Private */,
       9,    0,   61,    2, 0x08 /* Private */,

 // slots: parameters
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

void JQuadDominant2PureQuadsDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JQuadDominant2PureQuadsDialog *_t = static_cast<JQuadDominant2PureQuadsDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->enumFaces(); break;
        case 1: _t->getNewStrip(); break;
        case 2: _t->remeshStrip(); break;
        case 3: _t->getAllStrips(); break;
        case 4: _t->remeshAllStrips(); break;
        case 5: _t->refineAll(); break;
        case 6: _t->clearAll(); break;
        case 7: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JQuadDominant2PureQuadsDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JQuadDominant2PureQuadsDialog.data,
      qt_meta_data_JQuadDominant2PureQuadsDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JQuadDominant2PureQuadsDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JQuadDominant2PureQuadsDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JQuadDominant2PureQuadsDialog.stringdata0))
        return static_cast<void*>(const_cast< JQuadDominant2PureQuadsDialog*>(this));
    if (!strcmp(_clname, "Ui::QuadDominant2PureQuadsDialog"))
        return static_cast< Ui::QuadDominant2PureQuadsDialog*>(const_cast< JQuadDominant2PureQuadsDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JQuadDominant2PureQuadsDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 8;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
