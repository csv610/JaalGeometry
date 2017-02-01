/****************************************************************************
** Meta object code from reading C++ file 'EdgeAttributesDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "EdgeAttributesDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'EdgeAttributesDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JEdgeAttributesDialog_t {
    QByteArrayData data[13];
    char stringdata0[150];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JEdgeAttributesDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JEdgeAttributesDialog_t qt_meta_stringdata_JEdgeAttributesDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JEdgeAttributesDialog"
QT_MOC_LITERAL(1, 22, 8), // "setColor"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 8), // "setAlpha"
QT_MOC_LITERAL(4, 41, 8), // "setGlyph"
QT_MOC_LITERAL(5, 50, 12), // "setLineWidth"
QT_MOC_LITERAL(6, 63, 17), // "setCylinderRadius"
QT_MOC_LITERAL(7, 81, 16), // "setCylinderSides"
QT_MOC_LITERAL(8, 98, 12), // "checkDisplay"
QT_MOC_LITERAL(9, 111, 13), // "keyPressEvent"
QT_MOC_LITERAL(10, 125, 10), // "QKeyEvent*"
QT_MOC_LITERAL(11, 136, 1), // "e"
QT_MOC_LITERAL(12, 138, 11) // "closeDialog"

    },
    "JEdgeAttributesDialog\0setColor\0\0"
    "setAlpha\0setGlyph\0setLineWidth\0"
    "setCylinderRadius\0setCylinderSides\0"
    "checkDisplay\0keyPressEvent\0QKeyEvent*\0"
    "e\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JEdgeAttributesDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       9,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   59,    2, 0x08 /* Private */,
       3,    0,   60,    2, 0x08 /* Private */,
       4,    0,   61,    2, 0x08 /* Private */,
       5,    0,   62,    2, 0x08 /* Private */,
       6,    0,   63,    2, 0x08 /* Private */,
       7,    0,   64,    2, 0x08 /* Private */,
       8,    0,   65,    2, 0x08 /* Private */,
       9,    1,   66,    2, 0x08 /* Private */,
      12,    0,   69,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 10,   11,
    QMetaType::Void,

       0        // eod
};

void JEdgeAttributesDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JEdgeAttributesDialog *_t = static_cast<JEdgeAttributesDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->setColor(); break;
        case 1: _t->setAlpha(); break;
        case 2: _t->setGlyph(); break;
        case 3: _t->setLineWidth(); break;
        case 4: _t->setCylinderRadius(); break;
        case 5: _t->setCylinderSides(); break;
        case 6: _t->checkDisplay(); break;
        case 7: _t->keyPressEvent((*reinterpret_cast< QKeyEvent*(*)>(_a[1]))); break;
        case 8: _t->closeDialog(); break;
        default: ;
        }
    }
}

const QMetaObject JEdgeAttributesDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JEdgeAttributesDialog.data,
      qt_meta_data_JEdgeAttributesDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JEdgeAttributesDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JEdgeAttributesDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JEdgeAttributesDialog.stringdata0))
        return static_cast<void*>(const_cast< JEdgeAttributesDialog*>(this));
    if (!strcmp(_clname, "Ui::EdgeAttributesDialog"))
        return static_cast< Ui::EdgeAttributesDialog*>(const_cast< JEdgeAttributesDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JEdgeAttributesDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 9)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 9;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 9)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 9;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
